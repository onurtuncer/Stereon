# ADR-0001: Tolerance model

- **Status:** Proposed
- **Date:** 2026-09-28
- **Revised:** 2026-10-10 (coincidence is not identity; merge policy; three tolerance concepts)
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

### Rules

- **Units and size box.** Internal length unit is the millimetre. All geometry must lie inside a ±500 000 mm (±500 m) size box.
- **Linear resolution:** 1e-5 mm. **Angular resolution:** 1e-11 rad. Two points closer than the resolution cannot be distinguished by any *single* comparison. The resolution does **not** define an equivalence relation: coincidence tests answer one question about one pair; whether entities are merged is decided by the merge policy below.
- **Local tolerances.** Every `Edge` and `Vertex` stores a tolerance ≥ resolution. Native construction produces resolution-tight entities. A local tolerance may only grow through an explicit operation (import healing, sewing, boolean fallback) that logs the entity, the old and new value, and the reason, and reports the change in the `OperationResult` (ADR-0006).
- **Containment invariant.** A vertex tolerance sphere contains the ends of all its edges; an edge tolerance tube contains its 3D curve and all its pcurves mapped to 3D. The validity checker (ADR-0018) enforces this.
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
