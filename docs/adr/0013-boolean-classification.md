# ADR-0013: Boolean classification

- **Status:** Accepted
- **Date:** 2026-09-28
- **Revised:** 2026-10-10 (contact-type-aware cross-check; deterministic ray directions)
- **Revised:** 2026-10-14 (separate contact type from side transition; unresolved/fallback contract in ADR-0029)
- **Revised:** 2026-10-20 (the decisions on acceptance below: the seeding budget and the bounded winding-number fallback)
- **Accepted:** 2026-10-20 by Onur Tuncer, before the Boolean prototype that builds it for solids of planes and cylinders (`docs/planning/boolean-prototype-plan.md`, decision D7). The Phase 4 verification below is unchanged and remains the gate for the full domain.
- **Owner:** Onur Tuncer
- **Phase:** 4

## Context

After intersection and face splitting, a boolean must decide for every face fragment whether it lies inside, outside or on the boundary of the other solid. Misclassification is the second most common cause of wrong boolean results after intersection failures. The method must be robust for coincident faces, tangent contact and slightly imperfect input.

Not every intersection curve separates inside from outside. Where a face of one solid *crosses* the other solid's boundary, the fragments on either side have opposite states. Where it merely *touches* without crossing, both fragments are on the same side. Parallel normals alone do not prove this: higher-order tangential contact may still change sides (ADR-0029). A cross-check that demands opposite states everywhere rejects valid tangent models.

## Options

1. **Ray casting per fragment.** Simple; fails on rays grazing edges or faces, and is slow if done for every fragment.
2. **Generalised winding number.** Robust even for imperfect (non-watertight) input; expensive on curved B-rep surfaces.
3. **Hybrid propagation.** Classify one fragment per connected region with a robust point-in-solid test, then propagate the result across edges that were not intersected.

## Decision

Adopt option 3.

- **Contact and side-transition evidence.** Follow [ADR-0029](0029-intersection-and-classification-outcomes.md). Intersection output distinguishes transverse `Crossing`, `Tangential` contact and `Overlap` regions/boundaries. For each relevant adjacent-fragment pair, carry independently established SameSide, OppositeSide or Unresolved evidence. Tangential normals alone do not establish SameSide. Split branches where contact or transition changes; handle singularities and unresolved intervals explicitly.
- **Propagation.** After splitting, build fragment connectivity. Fragments across a non-intersection edge share classification only when the intersection set is complete for that neighborhood and no boundary transition is omitted. Across contact branches, propagate a same state only with SameSide evidence and toggle Inside/Outside only with OppositeSide evidence. Overlap boundaries use oriented overlap handling. Resolve or independently classify uncertain relations; never propagate through an unresolved transition as though it were SameSide.
- **Seeding.** For each connected region, classify one interior sample point with point-in-solid: ray casting with randomised directions, exact predicates for plane-dominant cases, and re-shooting when a ray hits within tolerance of an edge or vertex. Ray directions are a pure function of the input and the attempt number (ADR-0007): a fixed sequence, as the Boolean prototype's sixteen directions are, or a counter-based generator seeded with the region's lowest fragment index and the attempt number; no global or thread-local random state is used, so the same input gives the same rays regardless of scheduling.
- **Fallback.** If ray results disagree, a generalised winding number on a tessellation may be used only with the approximation/error and classification-insensitivity evidence required by ADR-0029. A fine mesh alone is insufficient. Refine or return an unresolved-classification failure when the requested assurance cannot be established; do not guess Outside or OnBoundary.
- **Coincident faces** require established trimmed-domain overlap and orientation evidence before classification as `OnSame` or `OnOpposite`, then are kept or removed according to the operation. Parallel normals or close supports alone do not establish that overlap (ADR-0029).
- **Cross-check.** Classification results are checked against established contact and side-transition evidence:

  | Boundary between two fragments of the same input face | Required relation |
  | --- | --- |
  | Established OppositeSide transition (including a regular transverse crossing) | opposite states (`Inside` / `Outside`) |
  | Established SameSide transition | same state |
  | `Overlap` boundary | one side `OnSame` or `OnOpposite`; the other side comes from its own justified classification |
  | not an intersection edge, with neighborhood completeness established | same state |
  | Unresolved contact or transition | refine/classify with adequate evidence, or fail; no forced relation |

  A violated established relation means intersection, splitting or seeding went wrong, and the operation fails with `KernelError::DegenerateResult` naming the edge and both fragments. Lack of evidence is a distinct unresolved-classification diagnostic, not a contradicted relation or a successful boundary classification.

### Decisions on acceptance (2026-10-20)

- **Seeding has a fixed budget.** Each region's seed point is classified by up to 16 rays, their directions from the counter-based generator above. A hit is decided with exact predicates: a ray meets a plane or a cylinder at the roots of a linear or quadratic equation, and whether the hit lies inside a trimmed face is a point-in-loop test against the face's line and circle pcurves. A ray that passes within tolerance of an edge or a vertex, or touches a surface tangentially, is discarded and the next direction tried.
- **Then the bounded winding number.** When no ray gives a clean answer, the seed point's generalised winding number decides, with an error bound: planar faces contribute their exact solid angles, and curved faces are tessellated with a known deviation bound that gives a bound on their contribution. The side is accepted only when the winding number lies further from ½ than its bound, so no error within the bound could change it (ADR-0029 §4). Otherwise the tessellation is refined, up to a budget.
- **Then `ErrorCode::Unresolved`**, naming the region, the rays tried and the last winding-number interval (ADR-0029 §6). The operation publishes nothing.

## Consequences

- Only a handful of point-in-solid queries per boolean; propagation does the rest.
- Coincident and tangent cases become explicit code paths instead of accidents, and tangent contact no longer trips the cross-check.
- Requires reliable fragment connectivity, which depends on correct face splitting, and reliable contact and side-transition evidence from ADR-0010/0029.
- Classification is reproducible across thread counts because the randomness is a pure function of the input.

## Verification

Phase 4 exit gate: ≥ 99.5% success on the 5 000-case boolean corpus; every failure returns a structured error; cross-check against OCCT and Manifold volumes within 1e-6 relative. The corpus includes a cylinder tangent to a plane face, internally tangent spheres, a cylinder tangent to another cylinder along a line, coplanar faces with `OnSame` and `OnOpposite` normals, and a tangent contact that becomes a crossing along the same branch; all must succeed, and the determinism run (ADR-0007) must give identical results.

The ADR-0029 contact corpus additionally tests ordinary tangency versus higher-order tangential crossing, ambiguity near a boundary, and tessellation fallbacks whose error prevents a reliable side decision. Unsupported/unresolved outcomes must be explicit until the declared supported domain includes those cases; they do not relax the mandatory cases above.
