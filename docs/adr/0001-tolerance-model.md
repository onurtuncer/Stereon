# ADR-0001: Tolerance model

- **Status:** Accepted
- **Date:** 2026-09-28
- **Revised:** 2026-10-10 (coincidence is not identity; merge policy; three tolerance concepts)
- **Revised:** 2026-10-14 (purpose-specific policy extension in ADR-0028)
- **Revised:** 2026-10-20 (resolutions and ceiling as `Context` values; the size box as the caller's obligation; the containment invariant checked by the Geometric validity rules; the Boolean and merge tests as conditions of the Boolean milestone)
- **Accepted:** 2026-10-20 by Onur Tuncer: the policy is built wherever code exists (local tolerances, explicit and reported growth, the ceiling, the epsilon lint, the containment check), and the tests that need Booleans are conditions of the Boolean prototype, whose failure reopens this ADR (see *Status of the verification*)
- **Owner:** Onur Tuncer
- **Phase:** 0

## Context

Every topological decision in a B-rep kernel (is this point on that surface? do these edges meet?) depends on a tolerance. Imported STEP data arrives with gaps of 1e-3 mm or more, while intersection curves are computed to 1e-7 mm. Kernels that answer this with scattered epsilons become unpredictable: the same boolean succeeds or fails depending on model scale. This decision shapes the topology store, every intersection routine and the boolean pipeline, so it must be fixed before any Phase 1 code.

A second trap is treating "within tolerance" as "the same entity". Proximity is not transitive: if A is within tolerance of B and B of C, A may still be far from C. A kernel that merges whatever is close produces a topology that depends on the order in which coincidences were discovered.

## Options

1. **Single global tolerance.** Simple, but imported models with larger gaps cannot be represented without healing everything to the global value.
2. **Exact geometry throughout** (rational arithmetic, exact curves). Removes tolerance issues in theory; impractical for NURBS surface intersections and far too slow.
3. **Global resolution plus local tolerances** (the Parasolid approach). A fixed model resolution governs all native geometry; edges and vertices may carry a larger local tolerance when data demands it.

## Decision

Adopt option 3.

### Three quantities that are kept apart

| Quantity | Meaning | Owner |
| --- | --- | --- |
| **Numerical error bound** | How far a computed value may be from the exact one. Attached to constructions (intersection points, fitted curves) as a certificate. | ADR-0002, ADR-0009, ADR-0010 |
| **Modelling tolerance** | The resolution and the local tolerances of vertices and edges: how far apart two entities may be and still be *treated* as coincident. | this ADR |
| **Healing budget** | How much a local tolerance may grow during one operation before it fails. The ceiling below, adjustable per `Context`. | this ADR, ADR-0006 |

A construction's error bound must be at most the modelling tolerance of the entity it produces; otherwise the entity's tolerance grows (logged) or the operation fails.

### Purpose-specific controls

[ADR-0028](0028-purpose-specific-numerical-tolerances.md) extends this policy with distinct relative, parametric, solver-convergence and approximation controls, explicit dimensional conversions and persistence/replay requirements. Numerical error bounds remain evidence, not adjustable epsilons. Its additions preserve this ADR's resolution, containment, tolerance-growth ceiling and merge policy.

### Rules

- **Units and size box.** Internal length unit is the millimetre. All geometry must lie inside a ±500 000 mm (±500 m) size box (`Context::GetSizeBox()`). The box is the caller's obligation: the kernel does not check it, and results for geometry outside it are unspecified, since the predicates' error bounds (ADR-0002) and the resolution assume it.
- **Linear resolution:** 1e-5 mm. **Angular resolution:** 1e-11 rad. Both, and the ceiling below, are `Context` values with these defaults (`WithResolution`, `WithAngularResolution`, `WithToleranceCeiling`), and a model's values are saved and restored with it, never defaulted on load (ADR-0028); the size box is fixed. Two points closer than the resolution cannot be distinguished by any *single* comparison. The resolution does **not** define an equivalence relation: coincidence tests answer one question about one pair; whether entities are merged is decided by the merge policy below.
- **Local tolerances.** Every `Edge` and `Vertex` stores a tolerance ≥ resolution. Native construction produces resolution-tight entities. A local tolerance may only grow through an explicit operation (import healing, sewing, boolean fallback) that logs the entity, the old and new value, and the reason, and reports the change in the `OperationResult` (ADR-0006).
- **Containment invariant.** A vertex tolerance sphere contains the ends of all its edges; an edge tolerance tube contains its 3D curve and all its pcurves mapped to 3D. The validity checker (ADR-0018) enforces this at its Geometric level: rule G01 for vertices and G02 for edges, which compare each pcurve on its surface (and the edge's 3D curve, where it has one) with the edge's reference trace at evenly spread parameters. `Commit` runs them when the `Context`'s validity level asks (the default) and always in checked builds. The check is sampled; certified distance bounds replace the samples when interval evaluation exists.
- **Ceiling.** Local tolerances above 0.01 mm make an operation fail with `KernelError::ToleranceExceeded` unless the caller raises the ceiling explicitly in the `Context`.
- **No hidden epsilons.** Topology and boolean code never compare with literal constants; they query `ToleranceOf(entity)` or `ctx.GetResolution()`.

### Merge policy

Whenever an operation may identify several vertices (or edges) as one — sewing, boolean vertex coincidence, import healing — it follows this policy, so the result depends only on the input and never on thread count or discovery order (ADR-0007):

1. **Candidates** are collected and processed in a fixed order: by input entity index, then by an operation-defined secondary key. Never by completion order.
2. **Representative.** The first candidate in that order opens a cluster and becomes its representative; the representative's position is kept, never averaged, so repeated operations do not drift.
3. **Joining.** A later candidate joins an existing cluster only if it lies within the pair tolerance of the **representative** (not of any member) and the cluster radius after joining stays within the healing budget. Otherwise it opens a new cluster, even if it is close to a non-representative member.
4. **Result.** The merged vertex takes the representative's position and a tolerance equal to the largest of the member tolerances and the cluster radius, clamped to ≥ resolution. Any growth is logged.
5. **Edges** merge by the same procedure, with the pair test being containment of one edge's curve in the other's tolerance tube at both ends and at the certified maximum-deviation points.

Single-linkage clustering ("merge anything that touches anything already merged") is not allowed anywhere in the kernel.

## Consequences

- Imported data can be represented faithfully and healed incrementally.
- Every intersection and classification routine must accept and propagate tolerances; entity pairs use the larger of their tolerances.
- Tolerance growth becomes observable and testable, which makes debugging failed booleans much easier.
- The merge policy makes the result of sewing depend on input ordering in the rare non-transitive cases (A–B close, B–C close, A–C not). That dependence is deterministic and documented; the alternative, an order-independent clustering, cannot bound the cluster radius.
- Very large (> 1 km) or very small (micro-scale) models are out of scope for v1.0.

## Verification

Prototype union/cut/common on boxes and cylinders by the end of Phase 0 (week 12), including a set of models with deliberately loosened edges (1e-4 to 1e-2 mm). Pass: all results valid, tolerance growth logged, no operation depends on a literal epsilon (enforced by `tools/check_epsilon.py` on `topo/` and `boolean/`).

Merge-policy test: sew faces whose corner vertices form chains at 0.8× the pair tolerance (A–B–C–D collinear). Pass: the clusters follow rules 1–4 exactly, no cluster radius exceeds the budget, and the result is byte-identical across thread counts and repeated runs; permuting the input order changes the result only as rule 1 predicts.

### Status of the verification

| Criterion | Status |
| --- | --- |
| Local tolerances grow only explicitly, and are reported | Met: `ShapeBuilder::SetTolerance(entity, value, reason)`; `Commit` reports every change in its `OperationResult` and warns `ToleranceGrown` (`tests/topo/OperationTests.cpp`) |
| The ceiling | Met: `Commit` fails with `ToleranceExceeded` above `Context::GetToleranceCeiling()` |
| No literal epsilons | Met: `tools/check_epsilon.py` on `topo/`, `boolean/` and `mesh/`, in lint and CTest |
| The containment invariant | Met, sampled: rules G01 and G02 (`tests/topo/GeometricValidityTests.cpp`) |
| Booleans on boxes and cylinders, with loosened edges | A condition of the Boolean prototype (the first Boolean milestone), not of Phase 0's other work: the prototype is accepted only when it passes, and a failure reopens this ADR |
| The merge-policy test | The same: a condition of the first operation that merges (sewing or Boolean vertex coincidence) |
