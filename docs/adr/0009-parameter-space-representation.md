# ADR-0009: Edge curves as surface-intersection curves

- **Status:** Proposed (to be accepted before Phase 1 topology work)
- **Date:** 2026-09-29
- **Revised:** 2026-10-10 (topological use separated from geometric kind; `Trim` and `Spatial` kinds; explicit primary side; certificates defined)
- **Owner:** Onur Tuncer
- **Phase:** 1
- **Supersedes:** the earlier draft of this ADR ("Parameter-space curves (pcurves)"), never accepted

## Context

A face is a region of a surface bounded by edges. Trimming, point classification, tessellation and face splitting in booleans work most reliably in each surface's (u, v) parameter space, so every edge needs a 2D curve (pcurve) on each face it bounds. The hard problem is keeping an edge's descriptions consistent: when an edge stores an independent 3D curve plus separate pcurves, the three can drift apart through approximation, healing and repeated operations, and tolerance problems follow.

Golovanov (*Geometric Modeling: The Mathematics of Shapes*, the design basis of the C3D kernel) describes a stronger representation: an edge between two faces **is** the intersection curve of their two surfaces, stored as the pair of surfaces plus a 2D curve on each, with the 3D curve derived from them. This ADR adopts that idea for the case it fits — two faces on two distinct surfaces — and names the cases it does not fit. We use the published concepts only; no C3D code or API is copied.

Two surfaces do not always determine a curve. Two faces split from the same plane share one surface, and the intersection of a surface with itself is the whole surface. A wire body has no face at all. Imported surfaces may not meet anywhere within the tolerance ceiling, so no refitting of pcurves can close the gap. These cases need their own representation rather than being forced into the intersection form.

## Options

1. **Independent 3D curve, pcurves computed on demand** by projection. Least storage; projection is slow and fails near singularities (cone apex, sphere poles).
2. **Independent 3D curve plus mandatory pcurves** (the earlier draft). Consistent at creation; three independent curves can drift apart later.
3. **Surface-intersection curve.** The edge's geometry is the pair (S₁, c₁), (S₂, c₂): two surfaces and a pcurve on each, sharing one parameter. The 3D curve is derived, with a cached 3D approximation for speed.
4. **Pcurves only, no 3D information.** Compact; no fast 3D queries, and the two pcurves' images have no stated agreement.

## Decision

Adopt option 3 as the representation for edges between two faces on distinct surfaces, with an exact-curve form for analytic cases and explicit kinds for the cases option 3 cannot express.

**1. Topological use and geometric kind are separate.** The number of coedges on an edge (0, 1, 2 or more, from the radial ring of ADR-0004) is a topological fact. The *curve kind* below says what defines the edge's geometry. Both are stored; the validity checker (ADR-0018) enforces the allowed combinations.

| Kind | What defines the geometry | Stored | Allowed coedge counts |
| --- | --- | --- | --- |
| `Intersection` | Two **distinct** surfaces and a pcurve on each | Two pcurves on a shared parameter range (one per coedge, more for non-manifold edges), the primary coedge, a cached 3D approximation with its certificate, tolerance | ≥ 2 |
| `Exact` | A closed-form 3D curve (plane ∩ plane → line, plane ∩ cylinder → line/circle/ellipse, coaxial cylinder ∩ plane → circle, …) | The exact 3D curve plus exact pcurves on every adjacent surface | ≥ 1 |
| `Trim` | A single pcurve on one surface shared by **both** adjacent faces (a face split for naming, meshing, colouring or import) | One pcurve; both coedges reference it with opposite sense; 3D curve S(c(t)) derived, cached | 2 |
| `Boundary` | One surface and one pcurve (sheet boundary, wire-bounded face, seam of a periodic surface) | One pcurve per coedge (a seam carries two on the same surface) | 1 (2 for a seam) |
| `Spatial` | An **independent 3D curve** is the truth; pcurves are projections of it onto the adjacent surfaces within tolerance | The 3D curve, optional pcurves per coedge, tolerance | ≥ 0 (0 for free wires) |
| `Degenerate` | A pole or apex (sphere pole, cone apex) | A pcurve with no 3D extent and the single 3D point | ≥ 1 |

`Spatial` exists for two reasons: wire bodies and sketch edges have no face, and imported edges whose surfaces do not meet within the ceiling cannot be represented as `Intersection` without lying about their tolerance. Healing may upgrade a `Spatial` edge to `Intersection` or `Exact` when the surfaces permit, and the upgrade is logged like any tolerance change.

**2. Shared parameterisation.** All pcurves of an edge use the same parameter t on the same range [t₀, t₁]. For an `Intersection` edge, S₁(c₁(t)) and S₂(c₂(t)) are points of the same edge at the same t. This makes the two sides directly comparable and removes reparameterisation from every downstream algorithm.

**3. Evaluation and the primary side.** The edge point of an `Intersection` edge at t is S₁(c₁(t)) from its **primary coedge**, whose ID is stored on the edge. The primary is chosen at creation (the coedge on the face with the lower input index in the creating operation, ADR-0007) and changes only through an explicit, logged operation. It is **not** derived from radial-ring order, so sewing or non-manifold edits that reorder the ring never change how an edge evaluates. The deviation |S₁(c₁(t)) − S₂(c₂(t))| is bounded by the edge tolerance (ADR-0001); the validity checker verifies it.

**4. Certificates.** A *certificate* is a bound that holds over the **whole** parameter range, not at sample points. Every `Intersection` and `Trim` edge carries one for its cached 3D curve, and every `Intersection` edge carries one for the gap between its sides. Certificates are computed by interval evaluation (ADR-0002) of |S₁(c₁(t)) − S₂(c₂(t))| and |approx(t) − S₁(c₁(t))| over each knot span of the pcurves, subdividing a span until the interval width is below the target. If the subdivision budget runs out before the target is met, the edge is **not** certified: the creating operation either widens the edge tolerance to the proven bound within the healing budget (logged, ADR-0001) or fails with `KernelError::CertificationFailed`. Sampling (for example 1 000 points per edge) is a regression check in tests, never a certificate.

**5. Cached 3D curve.** Each `Intersection` and `Trim` edge carries a B-spline approximation of its 3D curve with its certificate. It serves bounding boxes, BVH construction, display and quick distance queries. It is a cache, never the truth: it can be discarded and rebuilt from the surfaces and pcurves at any time, and the rebuild is deterministic (ADR-0007).

**6. Non-manifold edges.** An edge shared by three or more faces keeps one pcurve per coedge, all on the shared parameter; the primary pair defines the evaluation, and every other side must also lie within tolerance.

**7. Refinement.** Because the surfaces are exact, pcurves of an `Intersection` edge can be refitted at any time, without touching the faces, down to the tolerance **the surfaces permit**. Where the surfaces themselves do not meet (imported data), the residual gap is a lower bound on the achievable tolerance; the edge keeps that tolerance or is represented as `Spatial`. Tolerance growth remains logged as ADR-0001 requires.

**8. Construction.** Every edge-creating operation (surface–surface intersection per ADR-0010, sweep, boolean splitting, sewing, face splitting) outputs the edge kind, the pcurves on a shared parameter, the primary coedge where applicable and the certificates. Analytic intersections produce `Exact` edges whenever the closed form exists. Splitting a face along a curve in its own parameter space produces `Trim` edges.

**9. Data exchange.**
- **STEP export** writes `Exact` edges as their exact curves; `Intersection` edges as `INTERSECTION_CURVE` with master representation `PCURVE_S1` (the primary side) and both pcurves as associated geometry; `Trim` and `Boundary` edges as `SURFACE_CURVE` with master `PCURVE_S1`; `Spatial` edges as `SURFACE_CURVE` with master `CURVE_3D`, or as a bare curve for wires.
- **STEP import** of edges that carry only a 3D curve creates `Spatial` edges first; healing then computes pcurves by projection and converts each edge to `Intersection`, `Trim` or `Exact` where the surfaces permit, recording the resulting tolerance. Edges that cannot be converted within the ceiling stay `Spatial`.

## Consequences

- An edge's descriptions are consistent by construction: for `Intersection` edges the surfaces are the truth, the pcurves locate the edge on them, and the 3D curve is derived.
- Edge tolerance becomes a measured and certified quantity rather than a stored guess.
- Tessellation, trimming and classification work in 2D without projection.
- Evaluating an edge point costs a surface evaluation, which is slower than a plain curve; the cached 3D curve keeps hot queries fast, and cache invalidation must be handled wherever pcurves change.
- Six kinds instead of four, and a stored primary coedge: more cases in every edge-handling switch, in exchange for representing wires, same-surface splits and gappy imports honestly.
- ADR-0004's edge record stores the kind, the primary coedge, the cached 3D curve with its certificate and the tolerance; pcurves stay on the coedges.
- Every edge-creating algorithm, including surface–surface intersection, must produce synchronised pcurves and certificates, which raises the bar for ADR-0010.
- Interval evaluation of composed curves S(c(t)) is needed in `surface/` and `intersect/` before Phase 1 ends.

## Verification

- For all Phase 3 constructed solids and all healed STEP imports in the corpus: every `Intersection` edge carries a gap certificate ≤ its edge tolerance; every cached 3D curve carries a deviation certificate. Independently, a 1 000-sample regression check per edge never observes a deviation above the certified bound.
- Every analytic intersection in the corpus with a closed form produces an `Exact` edge; every same-surface face split produces a `Trim` edge; the wire corpus loads as `Spatial` edges with zero coedges and passes the validity checker.
- Discarding and rebuilding all cached 3D curves leaves every classification, mass property and tessellation unchanged (bit-identical results). Reversing the radial ring order of every edge in a test body changes no evaluated point.
- STEP round-trip (export, re-import) preserves edge kinds and keeps tolerances within 1.1× of the originals.
