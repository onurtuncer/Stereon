# ADR-0013: Boolean classification

- **Status:** Proposed (to be accepted before Phase 4)
- **Date:** 2026-09-28
- **Owner:** Onur Tuncer
- **Phase:** 4

## Context

After intersection and face splitting, a boolean must decide for every face fragment whether it lies inside, outside or on the boundary of the other solid. Misclassification is the second most common cause of wrong boolean results after intersection failures. The method must be robust for coincident faces, tangent contact and slightly imperfect input.

## Options

1. **Ray casting per fragment.** Simple; fails on rays grazing edges or faces, and is slow if done for every fragment.
2. **Generalised winding number.** Robust even for imperfect (non-watertight) input; expensive on curved B-rep surfaces.
3. **Hybrid propagation.** Classify one fragment per connected region with a robust point-in-solid test, then propagate the result across edges that were not intersected.

## Decision

Adopt option 3.

- After splitting, build fragment connectivity. Fragments connected across an edge that is not part of an intersection curve share classification.
- For each connected region, classify one interior sample point with point-in-solid: ray casting with randomised directions, exact predicates for plane-dominant cases, and re-shooting when a ray hits within tolerance of an edge or vertex.
- If ray results disagree, fall back to a generalised winding number computed on a fine tessellation.
- **Coincident faces** are detected during intersection and classified as `OnSame` or `OnOpposite` by comparing normals, then kept or removed according to the operation.
- Classification results are cross-checked: adjacent regions across an intersection curve must have opposite inside/outside states, or the operation fails with `KernelError::DegenerateResult`.

## Consequences

- Only a handful of point-in-solid queries per boolean; propagation does the rest.
- Coincident and tangent cases become explicit code paths instead of accidents.
- Requires reliable fragment connectivity, which depends on correct face splitting.

## Verification

Phase 4 exit gate: ≥ 99.5% success on the 5 000-case boolean corpus; every failure returns a structured error; cross-check against OCCT and Manifold volumes within 1e-6 relative.
