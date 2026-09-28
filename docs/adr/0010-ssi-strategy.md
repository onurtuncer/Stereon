# ADR-0010: Surface–surface intersection strategy

- **Status:** Proposed (to be accepted before Phase 2)
- **Date:** 2026-09-28
- **Owner:** Onur Tuncer
- **Phase:** 2

## Context

Surface–surface intersection (SSI) is the foundation of booleans, fillets and sewing, and the largest single source of kernel failures. It must find every intersection branch (including small closed loops), handle tangential and singular contact, and output 3D curves plus pcurves on both surfaces within a guaranteed tolerance.

## Options

1. **Lattice/sampling methods.** Simple; can miss small loops and gives no completeness guarantee.
2. **Algebraic methods (implicitisation).** Exact topology for low-degree surfaces; impractical for general NURBS.
3. **Subdivision to find branches, then marching.** Recursive subdivision with bounding volumes finds starting points on every branch; adaptive marching traces each branch.
4. **Hybrid:** closed-form analytic cases, then option 3 as the generic path.

## Decision

Adopt option 4.

- **Analytic table first.** Closed-form handlers for plane–plane, plane–quadric, quadric–quadric where known (e.g. coaxial cylinders, cylinder–sphere) and plane–torus; results are exact conics or lines where possible.
- **Generic path.** Subdivide both surfaces into Bézier patches with tight bounding volumes; prune with interval tests; use normal-cone tests to rule out small closed loops in a patch pair; collect seed points on every branch and on patch boundaries.
- **Marching.** Trace each branch with adaptive step size controlled by curvature and chordal error, refining each point onto both surfaces with Newton iteration.
- **Special points.** Detect tangential contact, branch points and singularities explicitly; they end or split branches instead of being stepped over.
- **Output.** Each branch becomes a B-spline 3D curve plus a pcurve on each surface (ADR-0009), with a certified maximum deviation.
- **Failure.** An unresolved case returns `KernelError::IntersectionFailed` with the patch pair involved; it is never silently skipped.

## Consequences

- Most real-world cases (planes, cylinders, cones) take fast, exact analytic paths.
- The generic path is expensive; BVH pruning and parallelism per patch pair (ADR-0007) are required.
- A dedicated SSI test corpus of ≥ 500 cases must be built before booleans start.

## Verification

Phase 2 exit gate: on the SSI corpus (tangent cylinders, near-coincident NURBS, small loops, pole contact), all branches found, no spurious branches, and a tolerance certificate on every curve.
