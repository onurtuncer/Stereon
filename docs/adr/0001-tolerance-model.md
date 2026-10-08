# ADR-0001: Tolerance model

- **Status:** Proposed
- **Date:** 2026-09-28
- **Owner:** Onur Tuncer
- **Phase:** 0

## Context

Every topological decision in a B-rep kernel (is this point on that surface? do these edges meet?) depends on a tolerance. Imported STEP data arrives with gaps of 1e-3 mm or more, while intersection curves are computed to 1e-7 mm. Kernels that answer this with scattered epsilons become unpredictable: the same boolean succeeds or fails depending on model scale. This decision shapes the topology store, every intersection routine and the boolean pipeline, so it must be fixed before any Phase 1 code.

## Options

1. **Single global tolerance.** Simple, but imported models with larger gaps cannot be represented without healing everything to the global value.
2. **Exact geometry throughout** (rational arithmetic, exact curves). Removes tolerance issues in theory; impractical for NURBS surface intersections and far too slow.
3. **Global resolution plus local tolerances** (the Parasolid approach). A fixed model resolution governs all native geometry; edges and vertices may carry a larger local tolerance when data demands it.

## Decision

Adopt option 3.

- **Units and size box.** Internal length unit is the millimetre. All geometry must lie inside a ±500 000 mm (±500 m) size box.
- **Linear resolution:** 1e-5 mm. Two points closer than this are the same point. **Angular resolution:** 1e-11 rad.
- **Local tolerances.** Every `Edge` and `Vertex` stores a tolerance ≥ resolution. Native construction produces resolution-tight entities. A local tolerance may only grow through an explicit operation (import healing, sewing, boolean fallback) that logs the entity, the old and new value, and the reason.
- **Containment invariant.** A vertex tolerance sphere contains the ends of all its edges; an edge tolerance tube contains its 3D curve and all its pcurves mapped to 3D. The validity checker (ADR-0004) enforces this.
- **Ceiling.** Local tolerances above 0.01 mm make an operation fail with `KernelError::ToleranceExceeded` unless the caller raises the ceiling explicitly in the `Context`.
- **No hidden epsilons.** Topology and boolean code never compare with literal constants; they query `ToleranceOf(entity)` or `ctx.GetResolution()`.

## Consequences

- Imported data can be represented faithfully and healed incrementally.
- Every intersection and classification routine must accept and propagate tolerances; entity pairs use the larger of their tolerances.
- Tolerance growth becomes observable and testable, which makes debugging failed booleans much easier.
- Very large (> 1 km) or very small (micro-scale) models are out of scope for v1.0.

## Verification

Prototype union/cut/common on boxes and cylinders by the end of Phase 0 (week 12), including a set of models with deliberately loosened edges (1e-4 to 1e-2 mm). Pass: all results valid, tolerance growth logged, no operation depends on a literal epsilon (enforced by a clang-tidy check on `topo/` and `boolean/`).
