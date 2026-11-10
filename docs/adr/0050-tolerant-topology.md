# ADR-0050: Tolerant topology

- **Status:** Proposed. On acceptance it supersedes two rules of ADR-0001: the pair rule (the larger of two tolerances becomes their sum) and the ceiling (an absolute ceiling becomes one relative to the model, capped by the absolute one). The rest of ADR-0001 stands.
- **Date:** 2026-11-09
- **Revised:** 2026-11-10 (numbered 0050; aligned with what ADR-0001, ADR-0009, ADR-0018 and the code already provide; the owner's decisions on the pair rule, the ceiling and pcurve authority)
- **Owner:** Onur Tuncer
- **Phase:** 2 (v0.2 "Intersections", interim steps 0.2-T1 and 0.2-T2), completed in v0.3 (0.3-T3)
- **Related:** ADR-0001, ADR-0004, ADR-0006, ADR-0009, ADR-0010, ADR-0018, ADR-0025, ADR-0027, ADR-0028, ADR-0029, ADR-0030, ADR-0044, ADR-0051

## Context

Boolean robustness in commercial kernels comes less from better surface/surface intersection (SSI) than from how they absorb the geometric gaps SSI leaves behind. Every floating-point intersection curve misses its neighbouring edges by some amount. A kernel that cannot represent that miss either fails the operation or silently produces a leaky shell.

The three incumbent approaches:

| Kernel | Mechanism | Outcome |
| --- | --- | --- |
| ACIS (since R4, late 1990s) | Explicit `TVERTEX` / `TEDGE` types with per-entity tolerance; pcurves are authoritative face boundaries; global `SPAresabs` for exact entities | Booleans succeed on imperfect geometry; healing APIs tighten tolerances afterwards |
| Parasolid (since mid-1990s) | Tolerant edges and vertices as per-entity data; strict rules on when tolerance may grow | Same outcome; widely regarded as the most robust Booleans in the industry |
| OpenCASCADE | Per-entity tolerance on every vertex, edge and face, but growth is unbounded and rarely tightened | Tolerances inflate through successive operations until coincidence tests become meaningless; a main source of OCCT Boolean failures |

Two lessons follow. First, per-entity tolerance is necessary: a single global epsilon cannot describe a model whose gaps vary by six orders of magnitude between a planar face and a fillet blend. Second, the concept alone does not set kernels apart. ACIS and OCCT both store per-entity tolerances; what differs is the discipline about how tolerance originates, propagates and is bounded. Retrofitting that discipline onto a topology designed without it is what OCCT has struggled with for twenty years.

Intersections are where tolerance originates: SSI fit error, a curve that misses an existing edge, a tolerant coincidence (ADR-0029 §6). v0.2 is building those services now, and Booleans (v0.4–v0.6) will consume them. The discipline has to be fixed before intersection results feed topology.

### What Stereon already has

Much of the field-level design exists, from ADR-0001 (accepted) and its implementation:

- `VertexRecord::Tolerance` and `EdgeRecord::Tolerance` (`src/topo/Records.cppm`); faces carry no tolerance.
- Growth only through `ShapeBuilder::SetTolerance(entity, value, reason)`, which records a `ToleranceChange {Old, New, Reason}` in the operation result (ADR-0030) and a `ToleranceGrown` diagnostic. The reason is free text today.
- Vertex clustering by ADR-0001's merge policy (`ShapeBuilder::MergeVertexCluster`).
- Edge kinds with their authoritative representation (ADR-0009, ADR-0025), certificates and cached 3D curves.
- The Geometric validity rules of ADR-0018 (`src/topo/ValidateGeometry.cpp`):
  - G01, a vertex tolerance sphere contains its edges' ends;
  - G02, an edge tolerance tube contains its 3D curve and its pcurves mapped to 3D;
  - G04, every tolerance lies between the resolution and the ceiling.

  The checks are sampled (16 intervals), not certified.
- An absolute ceiling, `Context` default 0.01 mm, enforced with `ToleranceExceeded`.

Missing:

- named growth causes;
- a ceiling relative to the model;
- tightening and healing;
- the sum rule for pair tests;
- a regression suite of tolerant models;
- certified, rather than sampled, containment checks.

## Options

| Option | Assessment |
| --- | --- |
| Single global epsilon (classic textbook B-rep) | Rejected. Cannot represent models whose gaps vary across six orders of magnitude; forces either a huge epsilon that merges real features or a tiny one that fails most Booleans on real data. |
| ACIS-style separate tolerant types (`TVERTEX`, `TEDGE`) | Rejected. Doubles the type surface every algorithm must handle; the distinction is a value, not a kind. ACIS did this for backward compatibility with an existing exact topology, a constraint Stereon does not have. |
| OCCT-style per-entity tolerance without growth rules | Rejected. Already shown in practice to decay into unusable models; the field alone does not buy robustness. |
| Exact or adaptive-precision arithmetic only (CGAL-style) | Rejected as the whole answer. It removes the gap problem for polyhedral and algebraic geometry, and Stereon uses it for decisions (ADR-0002). But NURBS/NURBS intersections have no closed form and still need approximation, and exact constructions throughout would miss the performance targets for interactive modelling. |
| Tolerance on faces as well as edges | Rejected. It makes the invariants harder to state and check. Every known production kernel expresses face gaps on edges, and nothing in the machining use case needs otherwise. |
| Defer until Booleans exist | Rejected. The cheapest moment to fix the rules is before the intersection services that create tolerances feed topology. |
| **Per-entity tolerance as a field, with named, bounded and reported growth (chosen)** | Builds on ADR-0001's fields and reporting; adds named causes, a model-relative ceiling, the sum rule and tightening. |

## Decision

Every vertex and edge carries its own tolerance, and tolerance growth is governed by explicit, auditable rules rather than left to individual algorithms.

1. **Tolerance is a field, not a type.** As ADR-0001 and ADR-0004 already provide, `Vertex` and `Edge` hold a tolerance, and an exact entity has the model resolution `res_abs` (the `Context` resolution, 1e-5 mm by default). One code path handles both cases; no algorithm branches on entity kind.
2. **Pcurves on every coedge; authority by edge kind.** Every coedge owns a pcurve on its face, planar faces included (`Spatial` edges of wire bodies, which have no face, excepted). Which representation is authoritative stays as ADR-0009 and ADR-0025 define it:
   - for `Intersection` and `Trim` edges, the pcurves are, and the 3D curve is a certified witness within the edge tolerance;
   - for `Exact` and `Spatial` edges, the 3D curve is.

   So a face stays watertight in its own domain even where 3D curves disagree, while exact analytic edges keep their exact curve.
3. **Tolerance grows only through a named cause.** `SetTolerance` remains the single write path. Its free-text reason is replaced by a cause from a closed set, with free text kept as detail:

   | Cause | Origin |
   | --- | --- |
   | `AbsorbGap` | an intersection curve not meeting an existing edge or vertex |
   | `MergeVertices` | ADR-0001's merge policy |
   | `PromoteEdge` | upgrading a `Spatial` edge to `Intersection` or `Exact` (ADR-0009 §1) |
   | `TolerantCoincidence` | ADR-0029 §6, the edge tolerances that cover a declared tolerant coincidence |
   | `CertificateWidening` | ADR-0009 §4 |
   | `HealingBudget` | ADR-0010, an uncertified branch widened within the healing budget |
   | `Import` | tolerances read from a file |

   No algorithm writes the field directly. Growth without a recorded cause is a bug, not a feature, and a CI check rejects any other writer.
4. **Tolerance has a ceiling relative to the model, capped absolutely.** Every entity tolerance is bounded by `res_max = min(1e-3 × the model's bounding diagonal, the Context's absolute ceiling)`. The absolute ceiling is ADR-0001's, 0.01 mm by default, so a large model never accepts millimetre gaps and a small one gets a tighter bound. An operation that would exceed `res_max` fails with `ToleranceExceeded` instead of producing a model that cannot be reasoned about. This supersedes ADR-0001's absolute-only ceiling.
5. **Pair tests use the sum.** Two entities coincide when their distance is at most `tol(a) + tol(b)`, never a bare epsilon and never `res_abs` alone. This supersedes ADR-0001's "entity pairs use the larger of their tolerances". The merge policy's pair tolerance (ADR-0001 step 3) and ADR-0029's tolerant coincidence follow it. One helper, `coincident(a, b)`, implements it, and code review rejects bare epsilons (ADR-0001's epsilon lint already enforces part of this).
6. **Tightening is a first-class operation.** `tighten(entity)` and `heal()` (below) lower tolerances to the achieved error. A successful Boolean leaves tolerances no worse than its inputs plus the SSI approximation error.
7. **Invariants are checked.** The existing Geometric rules G01, G02 and G04 are the three invariants below, with G04 reading `res_max`. Certified distance bounds replace their samples where interval evaluation exists (ADR-0001 already anticipates this). `validate()` is ADR-0018's checker at its Geometric level; no second checker is added.

### Data model

| Entity | Tolerance | Meaning |
| --- | --- | --- |
| `Vertex` | `Tolerance` (exists) | Radius of a ball around its point that contains the true meeting point of all incident edge curves |
| `Edge` | `Tolerance` (exists) | Radius of a tube around its reference trace that contains every adjacent coedge's pcurve mapped to 3D, and its 3D curve |
| `Coedge` | none | Owns its pcurve on its face |
| `Face` | none | Its surface is exact; its boundary is defined by its loops' pcurves |
| Model (`Context`) | `res_abs`, absolute ceiling (exist); `res_max` (new, derived) | Floor and ceiling for all entity tolerances |

### Invariants

$$\forall c \in \mathrm{coedges}(e):\quad \max_{t}\, \lVert S_{f(c)}(p_c(t)) - C_e(\phi_c(t)) \rVert \le \mathrm{tol}(e) \quad \text{(G02)}$$

$$\forall e \ni v:\quad \lVert C_e(t_{\mathrm{end}}) - P_v \rVert \le \mathrm{tol}(v) \quad \text{(G01)}$$

$$\mathrm{res\_abs} \le \mathrm{tol}(\cdot) \le \mathrm{res\_max} \quad \text{(G04)}$$

Here $S$ is the face surface, $p_c$ the pcurve, $C_e$ the edge's reference trace (its 3D curve or its primary pcurve on its surface, by edge kind) and $\phi_c$ the parameter map between pcurve and edge.

### Where tolerance originates

| Source | Typical magnitude | How it is absorbed |
| --- | --- | --- |
| SSI curve approximation (a fitted marched curve) | The fit error, at most the requested accuracy, by default `res_abs` | Edge tolerance set to the certified fit error at creation |
| Intersection curve not meeting an existing edge | Up to `res_max` | `AbsorbGap` on the shared vertex |
| Near-coincident surfaces in a Boolean | Variable | Declared a tolerant coincidence (ADR-0029 §6) within the pair tolerance; otherwise the edge tolerance is raised and flagged for healing |
| Imported geometry (STEP) with its own tolerances | Whatever the file says | Imported with `Import`; validated on import; entities above `res_max` rejected or healed |

### Healing and tightening

`tighten(entity)` refits what is not authoritative and lowers the tolerance to the achieved, certified error:

- for an `Intersection` or `Trim` edge, the cached 3D curve to its pcurves (ADR-0009 §7);
- for a vertex, its position, to the least-squares meeting point of its edges.

It never moves an authoritative representation, so it does not change the shape. Booleans call it on their output edges by default. `heal()` iterates until no tolerance can be lowered by more than a factor of two. Both are `ShapeBuilder` operators under ADR-0027, and both report their changes in the operation result.

## Consequences

**Positive**

- Booleans and SSI can succeed on imperfect inputs from their first implementation, instead of being rewritten once real models arrive.
- One code path for exact and tolerant entities; no algorithm branches on entity kind.
- Every tolerance in a model has a recorded, typed cause, so robustness regressions are diagnosable rather than mysterious.
- STEP import of models from ACIS, Parasolid and OCCT becomes a tolerance mapping, not a healing nightmare.
- The ceiling turns silent degradation into a loud, early failure, and scales with the model without allowing millimetre gaps.

**Negative**

- The sum rule is looser than ADR-0001's larger-of-two. The Boolean assembly and Euler split code, and the Boolean prototype plan's decision D2, follow it in a separate change. That change must rerun the Boolean prototype's fixtures and the merge-policy test.
- Every coincidence test is written against `tol(a) + tol(b)`; reviewers reject bare epsilons.
- Pcurves become mandatory on every face coedge, including planar faces where they are otherwise redundant: roughly one extra 2D curve per coedge.
- Certified containment checks cost more than today's sampled ones and stay debug or validation-level.
- Tightening after Booleans adds a refit pass whose cost scales with the number of new edges.

**Risks**

- Growth rules that are too strict cause failures on legitimate inputs; rules too loose reproduce OCCT's problem. `res_max` will need tuning against real models; the tolerant regression suite and the corpus are the way to tune it.
- The ACIS and Parasolid patents on tolerant modelling from the 1990s have expired, but any implementation detail taken from their documentation is reviewed before the SDK packaging ADR (ADR-0033) is finalised.

## Verification

Interim steps (roadmap revision 0.5):

| Step | Release | Work |
| --- | --- | --- |
| 0.2-T1 | v0.2 | The closed set of causes on `SetTolerance` and in `ToleranceChange`; `res_max` derived in `Context` and read by G04; the CI check that no code outside `ShapeBuilder::SetTolerance` writes a tolerance; `coincident(a, b)` with the sum rule, and the Boolean and Euler code moved to it with the prototype fixtures rerun |
| 0.2-T2 | v0.2 | Intersection results that create edges report the certified fit error and the gap that become tolerances, with their causes (ADR-0029 §6, ADR-0010) |
| 0.3-T3 | v0.3 | Pcurve synthesis for every face coedge; `tighten()` and `heal()`; certified G01/G02 bounds; the tolerant regression suite run through validation in CI |

Accepted when:

- [ ] The regression suite of tolerant models passes the Geometric checks: hand-built gaps at 1e-6, 1e-4 and 1e-2 of the bounding diagonal, the last rejected where it exceeds `res_max`.
- [ ] No tolerance is written outside `SetTolerance`, and every growth carries a cause from the closed set (a CI check).
- [ ] The Boolean prototype's fixtures and ADR-0001's merge-policy test pass under the sum rule.
- [ ] A Boolean between two bodies whose shared edge is displaced by up to `res_max` succeeds and leaves the output within tolerance (v0.4–v0.6 gate).
- [ ] The tolerances of Boolean outputs are no worse than the largest input tolerance plus the SSI fit error on the regression suite (the v0.6 tolerance-growth gate).
- [ ] A STEP file exported from OCCT with inflated tolerances imports, heals and validates (v0.8 gate).
