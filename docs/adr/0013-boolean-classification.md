# ADR-0013: Boolean classification

- **Status:** Proposed (to be accepted before Phase 4)
- **Date:** 2026-09-28
- **Revised:** 2026-10-10 (contact-type-aware cross-check; deterministic ray directions)
- **Owner:** Onur Tuncer
- **Phase:** 4

## Context

After intersection and face splitting, a boolean must decide for every face fragment whether it lies inside, outside or on the boundary of the other solid. Misclassification is the second most common cause of wrong boolean results after intersection failures. The method must be robust for coincident faces, tangent contact and slightly imperfect input.

Not every intersection curve separates inside from outside. Where a face of one solid *crosses* the other solid's boundary, the fragments on either side have opposite states. Where it merely *touches* (a cylinder resting on a plane, two spheres in internal tangency), both fragments are on the same side. A cross-check that demands opposite states everywhere rejects valid tangent models.

## Options

1. **Ray casting per fragment.** Simple; fails on rays grazing edges or faces, and is slow if done for every fragment.
2. **Generalised winding number.** Robust even for imperfect (non-watertight) input; expensive on curved B-rep surfaces.
3. **Hybrid propagation.** Classify one fragment per connected region with a robust point-in-solid test, then propagate the result across edges that were not intersected.

## Decision

Adopt option 3.

- **Contact types.** Every intersection edge produced for the boolean carries the contact type assigned by surface–surface intersection (ADR-0010): `Crossing` (the surfaces cross transversally along the curve), `Tangential` (normals parallel along the curve, surfaces on the same side) or `Overlap` (the curve bounds a region where the surfaces coincide). A branch whose contact type changes is split at the change point, so each edge has one type.
- **Propagation.** After splitting, build fragment connectivity. Fragments connected across an edge that is not an intersection edge share classification. Fragments connected across a `Tangential` edge also share classification.
- **Seeding.** For each connected region, classify one interior sample point with point-in-solid: ray casting with randomised directions, exact predicates for plane-dominant cases, and re-shooting when a ray hits within tolerance of an edge or vertex. Ray directions come from a counter-based generator seeded with the region's lowest fragment index and the attempt number (ADR-0007); no global or thread-local random state is used, so the same input gives the same rays regardless of scheduling.
- **Fallback.** If ray results disagree, fall back to a generalised winding number computed on a fine tessellation.
- **Coincident faces** are detected during intersection and classified as `OnSame` or `OnOpposite` by comparing normals, then kept or removed according to the operation.
- **Cross-check.** Classification results are checked against the contact types:

  | Edge between two fragments of the same input face | Required relation |
  | --- | --- |
  | `Crossing` | opposite states (`Inside` / `Outside`) |
  | `Tangential` | same state |
  | `Overlap` boundary | one side `OnSame` or `OnOpposite`; the other side's state comes from its own region seed and is not constrained |
  | not an intersection edge | same state |

  A violated relation means intersection, splitting or seeding went wrong, and the operation fails with `KernelError::DegenerateResult` naming the edge and both fragments.

## Consequences

- Only a handful of point-in-solid queries per boolean; propagation does the rest.
- Coincident and tangent cases become explicit code paths instead of accidents, and tangent contact no longer trips the cross-check.
- Requires reliable fragment connectivity, which depends on correct face splitting, and reliable contact typing from ADR-0010.
- Classification is reproducible across thread counts because the randomness is a pure function of the input.

## Verification

Phase 4 exit gate: ≥ 99.5% success on the 5 000-case boolean corpus; every failure returns a structured error; cross-check against OCCT and Manifold volumes within 1e-6 relative. The corpus includes a cylinder tangent to a plane face, internally tangent spheres, a cylinder tangent to another cylinder along a line, coplanar faces with `OnSame` and `OnOpposite` normals, and a tangent contact that becomes a crossing along the same branch; all must succeed, and the determinism run (ADR-0007) must give identical results.
