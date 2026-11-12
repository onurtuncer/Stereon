# ADR-0051: Face adjacency graph and canonical surface classification as kernel queries

- **Status:** Proposed
- **Date:** 2026-11-09
- **Revised:** 2026-11-10 (numbered 0051; the cache moved from the shell to an owned query workspace, after ADR-0005 and ADR-0036; related ADRs cited)
- **Owner:** Onur Tuncer
- **Phase:** 3 (v0.3 "Solids", interim steps 0.3-F1 and 0.3-F2), before the planar Booleans of v0.4, which consume the same queries
- **Related:** ADR-0005, ADR-0009, ADR-0013, ADR-0018, ADR-0025, ADR-0028, ADR-0029, ADR-0031, ADR-0036, ADR-0037, ADR-0038, ADR-0041, ADR-0044, ADR-0049, ADR-0050

Depends on ADR-0050 (tolerant topology): convexity classification reads edge tolerances. Complements ADR-0044 (tri-dexel stock) and ADR-0049 (machining collision and clearance queries). No feature recogniser is part of this ADR.

## Context

Four consumers inside and above Stereon need the same two questions answered about a B-rep. Which faces touch which, across which edges, and with what convexity? And what canonical surface is each face really?

| Consumer | What it needs | Why a bolt-on fails |
| --- | --- | --- |
| Boolean face splitting (v0.4) | Edge convexity, to decide which side of a split face survives; coincident-face detection | Convexity computed ad hoc inside the Boolean diverges from the recogniser's: two definitions of "tangent" |
| Rule-based feature recognition (AAG, hint-based) | A face adjacency graph with convexity labels; cylinder, plane and cone identity with axis and radius | An imported NURBS that is "really" a cylinder is invisible to every rule |
| Learned feature recognition (UV-Net, BRepNet) | The same graph as node and edge features, plus UV-sampled geometry | Rebuilding the graph outside the kernel re-implements adjacency, tolerance and orientation logic |
| Machining services (M7: ADR-0044 stock, ADR-0049 clearance) | Floor, wall and blend classification; feature volumes to decompose | Needs convexity consistent with the kernel's own notion, or the decomposition leaks |

In OpenCASCADE these queries exist only as scattered helpers (`BRepAdaptor`, `GeomConvert`, `ShapeAnalysis`), with no shared tolerance handling and no cache; every consumer rebuilds them. ACIS and Parasolid expose canonical surface identity and face–face relations as kernel API. That is one reason commercial feature recognisers (in NX and CAMWorks, for example) are built on them rather than on OCCT.

The edge-convexity question cannot be answered honestly without ADR-0050. At an edge with tolerance tol(e), the dihedral angle is known only to within an angular band that grows with tol(e). A kernel that returns "tangent" rather than "concave by 0.1 degrees" without reference to the edge tolerance is guessing.

Stereon today has the edge-to-face lookup (`FacesOf`) and traversal ranges (ADR-0036, ADR-0037), but no adjacency graph, dihedral or convexity code, and no recogniser from B-spline to analytic surfaces: `ToBSpline` converts only the other way.

## Options

| Option | Assessment |
| --- | --- |
| Leave adjacency and classification to the recogniser library | Rejected. Booleans need convexity before any recogniser exists; two implementations of "tangent" would drift apart; the recogniser would need tolerance internals it should not see. |
| Put feature recognition itself in the kernel | Rejected. Recognition rules are domain-specific (milling, turning, sheet metal) and evolve faster than the kernel; a learned recogniser has dependencies (an inference runtime) that do not belong in a geometry kernel. ADR-0049 keeps machining policy in applications for the same reason. |
| Classify convexity with a single angle, no band | Rejected. It produces confident wrong answers on tolerant edges: exactly the behaviour that makes OCCT-based recognisers fragile on imported data. |
| Store canonical type only when the surface is already analytic | Rejected. Imported STEP routinely delivers analytic shapes as NURBS; without fitting, every recogniser sees a freeform blob where there is a hole. |
| Compute queries on demand with no cache | Rejected. Booleans would rebuild the graph for every face pair; recognisers query the same faces repeatedly. |
| A cache stored inside the shape | Rejected. A frozen `Shape` is immutable (ADR-0005), and ADR-0036 keeps prepared adjacency out of the snapshot: "not stored in mutable Shape state". |
| Make the graph mutable by consumers | Rejected. Consumers would desynchronise it from the topology, and invalidation becomes impossible to reason about. |
| **Kernel queries over an immutable snapshot, results held in an owned, revision-keyed workspace (chosen)** | One definition of convexity and of canonical type; the snapshot stays immutable; invalidation follows from the snapshot's identity. |

## Decision

Stereon v0.3 ships two kernel-level queries, consumed by Booleans, recognisers and the machining services alike:

1. **Face adjacency graph (FAG).** For any shell, a graph whose nodes are faces and whose edges are the B-rep edges between them, each labelled with an assurance level and a convexity class drawn from `{convex, concave, tangent, knife, indeterminate}`. Convexity is computed from the dihedral angle along the edge and widened by the edge tolerance. An edge whose angular band straddles a class boundary is `indeterminate`, never silently rounded.
2. **Canonical surface classification.** For any face, the kernel reports whether its surface is, within a stated tolerance, a plane, cylinder, cone, sphere, torus or none of these. It distinguishes a sampled candidate from a whole-domain bounded classification. It also reports the canonical parameters (origin, axis, radius, half-angle). This is a property of the geometry, not of how the surface is stored: a NURBS imported from STEP that fits a cylinder within tol reports as a cylinder.

Both are kernel responsibilities, because both depend on tolerance semantics only the kernel owns, and because Booleans need them before any recogniser exists. Feature recognition itself (pocket, hole, slot, boss) stays out of the kernel, in a separate library above it.

Supporting queries ship alongside, because every known recogniser needs them:

- loop classification per face: the outer loop against the inner loops;
- a face role hint, `floor`, `wall`, `blend` or `other`, derived from the canonical type and the convexity of its boundary edges;
- vertex convexity: the sign pattern of the convexities of the edges meeting at a vertex;
- UV sampling of a face on a regular grid, returning positions, normals and the trim mask, for learned recognisers.

### Boundaries with existing decisions

- **Convexity is not classification evidence.** A dihedral class is not a side-transition proof. Boolean classification still propagates only on proved `SameSide`/`OppositeSide` evidence (ADR-0013, ADR-0029 §3; ADR-0028: "a normal-angle test alone does not establish side transition"). Booleans use convexity to choose which side of a split survives only where that evidence exists.
- **Recognition reports; it does not replace geometry.** `canonical()` returns a fit and its residual. It never swaps the stored surface for the analytic one. Replacing a surface is a separate, explicit healing operation under ADR-0027 and ADR-0050, with its tolerance consequences reported.
- **Not used by the validity checker.** ADR-0018's checker does not read adjacency caches, so the FAG is never validity evidence.
- **Geometric order is a query.** Radial coedge order is not angular (ADR-0037). The FAG's dihedral sampling is the geometric query that ADR-0037 says algorithms needing sector order must request. Non-manifold and `knife` edges use `ManifoldOppositeFace`'s checked preconditions, never an unqualified opposite face.
- **Face roles are query results, not attributes.** If an application stores them as attributes, ADR-0031's per-key policies apply.

### Query API

| Query | Signature (sketch) | Returns |
| --- | --- | --- |
| Adjacency graph | `Fag(const ShapeQueries&, ShellId) -> const FaceGraph&` | Nodes are face handles; edges are B-rep edge handles with a `Convexity` label, assurance and the angular band `[theta_min, theta_max]` |
| Canonical surface | `Canonical(const ShapeQueries&, FaceId, double tol) -> CanonicalSurface` | A variant of `Plane / Cylinder / Cone / Sphere / Torus / Freeform`, each with its parameters, assurance (`SampledCandidate`, `Bounded`, `Unresolved`), covered trim domain and sampled residual or certified whole-domain deviation bound |
| Loop classes | `Loops(const ShapeQueries&, FaceId) -> span<const LoopInfo>` | The outer/inner flag, area sign and bounding box in UV |
| Face role | `Role(const ShapeQueries&, FaceId) -> FaceRole` | `floor / wall / blend / other`, with the convexity evidence used |
| Vertex convexity | `VertexConvexity(const ShapeQueries&, VertexId) -> VertexClass` | The sign pattern of the incident edges' convexities |
| UV sample | `SampleUv(const Shape&, FaceId, int nu, int nv) -> UvGrid` | Positions, normals and trim mask, row-major |

`ShapeQueries` is the owned workspace below. No query mutates topology. Names follow ADR-0017.

### Convexity classification

Take edge $e$, shared by faces $A$ and $B$, at parameter $t$ along it. Let $n_A$ and $n_B$ be the outward normals, and $d$ the edge tangent oriented by face $A$'s coedge. Then

$$\theta(t) = \operatorname{atan2}\big( (n_A \times n_B)\cdot d,\; n_A \cdot n_B \big)$$

Sampled at $N$ points along the edge, this gives only a sampled candidate interval $[\theta_{\min}, \theta_{\max}]$. The band is then widened by the angular uncertainty the edge tolerance induces. Let $\mathrm{tol}(e)$ be the edge tolerance (ADR-0050) and $\kappa$ the larger of the two faces' principal curvatures normal to the edge. Each face's normal direction is then uncertain by roughly $\kappa \cdot \mathrm{tol}(e)$ radians, so the band grows by $2\kappa\,\mathrm{tol}(e)$ on each side. An edge whose widened band is wider than a class makes the edge `indeterminate`.

| Class | Condition on the widened band |
| --- | --- |
| `convex` | entirely in $(\varepsilon_{\tan},\ \pi - \varepsilon_{\text{knife}})$ |
| `concave` | entirely in $(-\pi + \varepsilon_{\text{knife}},\ -\varepsilon_{\tan})$ |
| `tangent` | entirely within $\pm\varepsilon_{\tan}$ of 0 |
| `knife` | within $\varepsilon_{\text{knife}}$ of $\pm\pi$ |
| `indeterminate` | the band crosses a boundary |

Both $\varepsilon_{\tan}$ and $\varepsilon_{\text{knife}}$ default to 1 degree. They are named, purpose-specific angular allowances in the `Context`, under ADR-0028, and are persisted with it. A Boolean requiring a definite decision must obtain independent side-transition evidence or return an explicit unresolved/failure outcome. There is no `ConvexityOr` escape hatch that upgrades a heuristic to evidence. A recogniser treats `indeterminate` as a signal to look at neighbouring edges.

Sampling and the first-order estimate above are candidate generators only.
A bounded class requires an enclosure of the dihedral angle over the entire
edge interval, including normal/tangent regularity and angular branch cuts.
Curvature bounds must cover the relevant neighborhood; a sampled curvature
or the approximate `2*kappa*tol` formula is not a certificate. If this cannot
be established, report `SampledCandidate` or `Unresolved`, never a bounded
whole-edge class. Side-transition evidence remains separately required.

### Canonical surface recognition

For a face whose stored surface is already analytic, classification is a type check and a parameter copy. For a freeform surface:

1. Sample the surface on a UV grid inside the trim.
2. Fit, in order of simplicity:
   - a plane, by PCA of the samples;
   - a sphere, by an algebraic sphere fit;
   - a cylinder: the axis from PCA of the normal field, the radius by least squares;
   - a cone: the apex and half-angle from the normals;
   - a torus: two radii about the fitted axis.
3. Report the first fit whose maximum sampled residual is at most `tol` as `SampledCandidate`. Failure to find a fit is not proof of non-canonicity; report the search limitations.
4. Verify the accepted fit against the face's boundary curves, not only interior samples, so that a cylinder patch that bulges at its edge is not accepted.

Return sampled residuals separately from certified bounds. Boundary checks
alone do not exclude an interior bulge. Promotion to `Bounded` requires a
whole-trim-domain deviation bound, with method, domain, numerical policy and
limitations recorded. The default `tol` is the model's `res_abs`.

Analytic dispatch on a fitted candidate, even when explicitly requested, is
only a candidate-generation optimization. A correctness-critical result must
be verified against the original geometry, including completeness, trimming,
contact classification and accuracy; a small position-deviation bound alone
does not preserve intersection topology near contact. Without that evidence,
use the original-geometry solver or return an explicit unresolved outcome.
Recognition never silently changes the stored surface.

### Caching and invalidation

Results live in a `ShapeQueries` workspace that the caller owns, as ADR-0036 does for prepared adjacency (`TopologyQueryPreparation`) and ADR-0049 for acceleration data ("an explicitly owned workspace keyed by shape revision"). The snapshot itself stays immutable (ADR-0005).

- A workspace is bound to one frozen snapshot, identified by its store identity, and to the numerical context of its queries: the tolerances, the angular allowances and the requested fit tolerance.
- Results are computed lazily and kept for the workspace's lifetime.
- A different snapshot (any edit publishes a new one, ADR-0040) or a different numerical context means a different workspace, so stale results cannot be observed (ADR-0041: a relevant change "invalidates dependent evidence").
- Booleans build a workspace for their inputs and one for their output.
- A builder-side query during editing uses ADR-0038's topology epoch for connectivity. Results that depend on geometry also need a geometry epoch, which ADR-0038 does not define yet. 0.3-F1 adds it, or limits builder-side queries to topology-only results.

## Consequences

**Positive**

- One definition of convexity and one of canonical type across Booleans, recognisers and machining services.
- Recognisers, rule-based or learned, become thin libraries over a stable API, developed, versioned and distributed independently of the kernel.
- Imported STEP geometry behaves like native analytic geometry for every downstream consumer.
- `indeterminate` makes tolerance-induced ambiguity visible instead of hiding it in a wrong label.
- Immutable snapshots stay free of mutable cache state; thread safety follows from ADR-0005 and ADR-0007.

**Negative**

- Canonical fitting on every freeform face costs sampling and a least-squares solve per face at first query. The workspace amortises it, but it is noticeable on first load of large assemblies.
- A workspace adds memory roughly proportional to the edge count plus one small record per face.
- Callers must keep the workspace with the snapshot. A workspace kept for a superseded snapshot answers for that snapshot, correctly, but not for the new one.

**Risks**

- Class widths that suit prismatic parts may misclassify turned or blended parts; the defaults will need tuning against the corpus.
- Canonical fitting can accept a near-cylinder that a downstream CAM operation then treats as exactly cylindrical. The returned residual exists so that consumers can refuse; the risk is consumers ignoring it.
- The interface between kernel and recogniser library is a stability commitment: changing `FaceGraph` after v1.0 breaks every recogniser built on it.

## Verification

Interim steps (roadmap revision 0.5):

| Step | Release | Work |
| --- | --- | --- |
| 0.3-F1 | v0.3 | The `ShapeQueries` workspace keyed by snapshot and context; dihedral sampling and the widened-band classifier, with `Convexity`, the band and `ConvexityOr`; `FaceGraph` and `Fag`; loops and vertex convexity |
| 0.3-F2 | v0.3 | `Canonical` for analytic surfaces (a type check), then the freeform fitting chain (plane, sphere, cylinder, cone, torus) with boundary verification; `Role`; `SampleUv`; a minimal reference recogniser (through-hole, blind hole, rectangular pocket) in a separate library, the API's first consumer and proof that no kernel internals leak |
| 0.4 gate | v0.4 | Planar Boolean face splitting reads F1's convexity, together with ADR-0013 evidence |

Accepted when:

- [ ] Every edge in the synthetic corpus is classified as its construction intended, or `indeterminate` where an injected gap makes it so.
- [ ] A STEP export from OCCT of a cylinder as a NURBS surface reports `Cylinder`, with axis and radius within `res_abs` (with the v0.8 exchange gate; before it, the same surface built in-kernel).
- [ ] No stale result is observable after any Boolean or healing operation, verified by a fuzz test that edits shapes, publishes new snapshots and queries both the old and the new workspace.
- [ ] The reference recogniser finds every hole and pocket in the MFCAD subset used for testing, through the public query API only.
- [ ] On a 10 000-face shell, the first query takes under 1 s and repeated queries under 1 ms, measured as ADR-0043 specifies.

### Assurance regressions required before acceptance

- A NURBS patch with an interior bulge between grid samples and an unchanged
  boundary must not receive a bounded plane/cylinder classification.
- An edge whose dihedral changes class between samples must not receive a
  bounded single class; include angle wrapping and singular normals.
- Candidate analytic intersections must be checked against the original
  surface, including a near-tangent case with different component counts.
