# ADR-0010: Surface–surface intersection strategy

- **Status:** Accepted in part: the analytic table (2026-10-20); the generic path stays Proposed until Phase 2
- **Date:** 2026-09-28
- **Revised:** 2026-10-10 (completeness and accuracy as separate obligations; contact types)
- **Revised:** 2026-10-14 (dimension-aware outcomes and independent side-transition evidence in ADR-0029)
- **Revised:** 2026-10-20 (the decisions on acceptance below: the analytic table accepted, a cylinder–cylinder handler for skew and intersecting axes, a Tier A corpus for the Boolean prototype)
- **Accepted (analytic table):** 2026-10-20 by Onur Tuncer, before the Boolean prototype that builds on it (`docs/planning/boolean-prototype-plan.md`, decision D7). The generic path (subdivision, marching, normal-cone pruning) and the Phase 2 verification below stay Proposed until Phase 2.
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
- **Two separate obligations.** SSI must discharge both, and a failure of either is reported as such:
  1. **Completeness** — every branch is found. A patch pair is discarded only when a bounding-volume or interval test *proves* separation, or when the normal-cone test proves at most one branch and that branch has been found. A patch pair that cannot be resolved within the subdivision budget returns `KernelError::IntersectionFailed` with the patch pair; it is never silently skipped.
  2. **Accuracy** — every output curve carries a certificate (ADR-0009) bounding its deviation from both surfaces over the whole parameter range, computed by interval evaluation, not sampling. A branch that cannot be certified to the requested tolerance returns `KernelError::CertificationFailed` or is widened within the healing budget (logged).
- **Contact and side transition.** Follow [ADR-0029](0029-intersection-and-classification-outcomes.md): distinguish transverse `Crossing` branches, `Tangential` branches and `Overlap` regions/boundaries. Normal bounds help establish regular contact type, but parallel normals alone prove neither coincidence nor a same-side transition. Carry separate SameSide/OppositeSide/Unresolved evidence for the relevant adjacent fragments. Split where contact or transition changes; unresolved domains require refinement or explicit failure before downstream use. Overlap regions carry their extent and parameter correspondence, with boundaries where present.
- **Output.** A complete intersection may contain isolated points, curves and overlap regions, or be provably empty over the requested domain (ADR-0029). Curves use exact analytic geometry where available, otherwise the ADR-0009 surface/pcurve representation and certified 3D approximation; include parameter domains, certificates and contact/transition evidence. Coincident supports do not imply overlapping trims. Partial components are diagnostic evidence on failure, never a complete successful result.
- **Inputs.** Both surfaces must be tier 2 (`BoundedGeometry`, ADR-0003) for the generic path; tier 1 surfaces return `KernelError::NotSupported`.

### Decisions on acceptance (2026-10-20)

- **The analytic table is accepted; the generic path is not yet.** The Boolean prototype builds the handlers for plane–plane, plane–cylinder and cylinder–cylinder. Subdivision, marching and normal-cone pruning stay Proposed, and are decided with the Phase 2 work that builds them.
- **Cylinder–cylinder in full, including skew and intersecting axes.** Parallel axes have closed forms (lines, or a coaxial overlap). Skew and intersecting axes give a quartic with no closed form, which RS-0001 leaves as Tier B; the prototype handles it with a dedicated handler rather than `NotSupported`. On one cylinder, at angle u, the other cylinder's equation is a quadratic in the height v, so the curve is the two branches v±(u) where its discriminant is not negative, meeting at the branch points where it vanishes, found exactly. Each branch's pcurves on both cylinders are B-splines fitted to the closed form and certified by interval evaluation of the residual (decision D3), and the edges are ADR-0009 `Intersection` edges with a cached 3D curve. Tangency between the cylinders (a double root along a branch) is decided exactly and reported as such.
- **The corpus.** The prototype needs a complete Tier A analytic corpus, plus the cylinder–cylinder cases above, with analytic ground truth and OCCT agreement, before its handlers merge: every plane–plane, plane–cylinder and cylinder–cylinder configuration, disjoint, tangent, crossing and coincident, near-degenerate on either side of each decision, and loosened. The 500-case corpus stays the gate of the generic path.

## Consequences

- Most real-world cases (planes, cylinders, cones) take fast, exact analytic paths.
- The generic path is expensive; BVH pruning and parallelism per patch pair (ADR-0007) are required.
- A dedicated SSI test corpus of ≥ 500 cases must be built before booleans start.

## Verification

Phase 2 exit gate: on the SSI corpus (tangent cylinders, near-coincident NURBS, small loops, pole contact), all branches found, no spurious branches, a certificate on every curve and correct contact/side-transition evidence on every relevant branch. Include isolated point and full/partial overlap outcomes from ADR-0029. Completeness is checked against ground truth where it exists (analytic cases) and against OCCT and a dense-sampling oracle elsewhere; accuracy is checked by a 1 000-sample regression run per curve, which must never exceed the certified bound. Cases that exhaust the budget must fail with the right error code rather than return a partial result.
