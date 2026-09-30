# ADR-0009: Edge curves as surface-intersection curves

- **Status:** Proposed (to be accepted before Phase 1 topology work)
- **Date:** 2026-09-29
- **Owner:** Onur Tuncer
- **Phase:** 1
- **Supersedes:** the earlier draft of this ADR ("Parameter-space curves (pcurves)"), never accepted

## Context

A face is a region of a surface bounded by edges. Trimming, point classification, tessellation and face splitting in booleans work most reliably in each surface's (u, v) parameter space, so every edge needs a 2D curve (pcurve) on each face it bounds. The hard problem is keeping an edge's descriptions consistent: when an edge stores an independent 3D curve plus separate pcurves, the three can drift apart through approximation, healing and repeated operations, and tolerance problems follow.

Golovanov (*Geometric Modeling: The Mathematics of Shapes*, the design basis of the C3D kernel) describes a stronger representation: an edge between two faces **is** the intersection curve of their two surfaces, stored as the pair of surfaces plus a 2D curve on each, with the 3D curve derived from them. This ADR adopts that idea. We use the published concepts only; no C3D code or API is copied.

## Options

1. **Independent 3D curve, pcurves computed on demand** by projection. Least storage; projection is slow and fails near singularities (cone apex, sphere poles).
2. **Independent 3D curve plus mandatory pcurves** (the earlier draft). Consistent at creation; three independent curves can drift apart later.
3. **Surface-intersection curve.** The edge's geometry is the pair (S₁, c₁), (S₂, c₂): two surfaces and a pcurve on each, sharing one parameter. The 3D curve is derived, with a cached 3D approximation for speed.
4. **Pcurves only, no 3D information.** Compact; no fast 3D queries, and the two pcurves' images have no stated agreement.

## Decision

Adopt option 3 as the general edge representation, with an exact-curve form for analytic cases.

**1. Edge curve kinds.** Every edge has one of four kinds:

| Kind | Used for | Stored |
| --- | --- | --- |
| `Intersection` | An edge between two faces (the general case) | Two surface references, two pcurves on a shared parameter range, a cached 3D approximation, tolerance |
| `Exact` | Both sides analytic and the intersection has a closed form (plane ∩ plane → line, plane ∩ cylinder → line/circle/ellipse, coaxial cylinder ∩ plane → circle, …) | The exact 3D curve plus exact pcurves on both surfaces |
| `Boundary` | An edge with one face (sheet boundary, wire-bounded face), or a seam edge on a periodic surface | One surface and one pcurve per coedge (a seam carries two pcurves on the same surface) |
| `Degenerate` | A pole or apex (sphere pole, cone apex) | A pcurve with no 3D extent and the single 3D point |

**2. Shared parameterisation.** Both pcurves of an `Intersection` edge use the same parameter t on the same range [t₀, t₁]. S₁(c₁(t)) and S₂(c₂(t)) are points of the same edge at the same t. This makes the two sides directly comparable and removes reparameterisation from every downstream algorithm.

**3. Evaluation.** The edge point at t is S₁(c₁(t)) from the primary side (the first coedge of the radial ring, ADR-0004). The deviation |S₁(c₁(t)) − S₂(c₂(t))| is bounded by the edge tolerance (ADR-0001); the validity checker verifies it.

**4. Cached 3D curve.** Each `Intersection` edge carries a B-spline approximation of its 3D curve with a certified maximum deviation. It serves bounding boxes, BVH construction, display and quick distance queries. It is a cache, never the truth: it can be discarded and rebuilt from the surfaces and pcurves at any time.

**5. Non-manifold edges.** An edge shared by three or more faces keeps one pcurve per coedge, all on the shared parameter; the primary pair defines the evaluation, and every other side must also lie within tolerance.

**6. Refinement.** Because the surfaces are exact, pcurves can be refitted to a tighter tolerance at any time (for example when healing or before a boolean) without touching the faces. Tolerance growth remains logged as ADR-0001 requires.

**7. Construction.** Every edge-creating operation (surface–surface intersection per ADR-0010, sweep, boolean splitting, sewing) outputs the edge kind, both pcurves on a shared parameter and the tolerance certificate. Analytic intersections produce `Exact` edges whenever the closed form exists.

**8. Data exchange.**
- **STEP export** writes `Exact` edges as their exact curves; `Intersection` edges as STEP `INTERSECTION_CURVE` or `SURFACE_CURVE` entities with the 3D approximation and both pcurves as associated geometry.
- **STEP import** of edges that carry only a 3D curve computes pcurves by projection and converts each edge to `Intersection` (or `Exact`, when recognised) during healing, recording the resulting tolerance.

## Consequences

- An edge's descriptions are consistent by construction: the surfaces are the truth, the pcurves locate the edge on them, and the 3D curve is derived.
- Edge tolerance becomes a measured quantity (the gap between the two sides) rather than a stored guess.
- Tessellation, trimming and classification work in 2D without projection.
- Evaluating an edge point costs a surface evaluation, which is slower than a plain curve; the cached 3D curve keeps hot queries fast, and cache invalidation must be handled wherever pcurves change.
- ADR-0004's edge record changes: it stores the curve kind, the cached 3D curve and the tolerance; pcurves stay on the coedges.
- Every edge-creating algorithm, including surface–surface intersection, must produce synchronised pcurves, which raises the bar for ADR-0010.
- Importers need a projection-and-conversion step for foreign data.

## Verification

- For all Phase 3 constructed solids and all healed STEP imports in the corpus: every `Intersection` edge satisfies |S₁(c₁(t)) − S₂(c₂(t))| ≤ edge tolerance at 1 000 samples per edge, and the cached 3D curve stays within its certified deviation.
- Every analytic intersection in the corpus with a closed form produces an `Exact` edge.
- Discarding and rebuilding all cached 3D curves leaves every classification, mass property and tessellation unchanged (bit-identical results).
- STEP round-trip (export, re-import) preserves edge kinds and keeps tolerances within 1.1× of the originals.
