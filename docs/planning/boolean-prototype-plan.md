# Boolean prototype plan (ADR-0001, the last Phase 0 gate)

**Owner:** Onur Tuncer
**Opened:** 2026-10-20
**Status:** decided 2026-10-20 (section 4); the work starts with M0, the review of the ADRs it builds on.

The Phase 0 exit gate has one open line: ADR-0007 accepted, which waits for a Boolean prototype that can test its speed-up and cross-platform checks. ADR-0001 is accepted with its Boolean and merge tests as conditions of the same prototype. This plan sets out what the prototype is, what it must prove, and the order of the work.

## 1. What the prototype is

Union, cut and common of solids bounded by planes and cylinders: boxes, cylinders and the solids they make together. Cylinder pairs with skew or intersecting axes (RS-0001 Tier B) are in scope too, through a dedicated handler (ADR-0010's acceptance). Out of scope: other surfaces, sewing of imported data and blends.

It must pass:

| Check | Source |
| --- | --- |
| Union, cut and common on boxes and cylinders, every result valid (Structural and Geometric), every tolerance growth logged, no literal epsilon in `topo/` or `boolean/` | ADR-0001 verification |
| The same with deliberately loosened inputs, 1e-4 to 1e-2 mm | ADR-0001 verification |
| The merge-policy chain test at 0.8× the pair tolerance, byte-identical across thread counts and runs | ADR-0001 verification |
| Every Tier A case of `corpus/analytic-ssi/` decided as RS-0001 requires | TODO.md §10, RS-0001 |
| Identical output at 1, 2, 8 and 32 threads; a speed-up of at least 5× on 8 cores on a many-face fixture; the cross-platform (G2) comparison | ADR-0007 verification |
| Splitting a face and an edge within ADR-0005's copy budget | ADR-0005 verification |
| Volumes, areas and face counts against OCCT's `BRepAlgoAPI` on shared cases | TODO.md §9 |

## 2. The pipeline

1. **Pair finding.** Candidate face pairs by bounding box, in input order.
2. **Surface intersection.** Analytic handlers for plane–plane, plane–cylinder and parallel cylinder–cylinder, returning ADR-0029 outcomes: proved empty, isolated points, crossing or tangential branches, overlap regions, full coincidence. Decisions are signs of polynomials in the input doubles through `stereon::robust` (RS-0001 G5).
3. **Face clipping.** Each branch is restricted to both trimmed faces: curve against loop in each face's parameter space.
4. **Imprinting.** Edges split where branches cross them, then the intersection edges are added and faces split along them. These are the controlled operators of ADR-0027, inside a transaction: a builder copy whose pages are shared until written, discarded on failure.
5. **Classification.** Fragments are grouped into connected regions. One seed per region is classified by a deterministic ray cast (ADR-0013), then the state is propagated across edges with SameSide or OppositeSide evidence (ADR-0029 §3). An unresolved transition is never crossed, and is never guessed.
6. **Assembly.** Fragments are kept per operation (OnSame and OnOpposite by the usual rules). Vertices and edges are merged by ADR-0001's representative policy, lower-dimensional contacts are dropped, and the result is published through `Commit` with a two-input History (ADR-0030).

## 3. Milestones

Each milestone merges on its own, with tests. Estimates are working days.

| # | Milestone | Content | Est. |
| --- | --- | --- | --- |
| M0 | ADR review | ADRs 0010 (SSI), 0013 (classification), 0027 (controlled operators), 0029 (outcomes) and 0030 (lineage) reviewed with the owner and accepted, each narrowed to the prototype's scope, before the milestones they govern. | 2 |
| M1 | Primitives | `Stereon.Build`: `MakeBox` and `MakeCylinder` solids with `Exact` edges and exact pcurves, valid at the Geometric level; the test bodies' builders reused. | 2 |
| M2 | Outcome types | `SurfaceSurfaceIntersection` grows into ADR-0029's outcome: kind, parameter domains on both surfaces, contact, side-transition evidence, a proved-empty flag, overlap correspondence, and the snaps applied. | 2 |
| M3 | Analytic handlers | Plane–plane; plane–cylinder (disjoint, tangent line, two lines, circle, ellipse); parallel cylinder–cylinder (disjoint, tangent, two lines, coincident). Exact decisions, `Exact` curves and pcurves on one parameter. Runs on `corpus/analytic-ssi/` and against OCCT. | 6 |
| M3b | Cylinder–cylinder, skew and intersecting axes | The quartic as two branches v±(u) on one cylinder, branch points found exactly, pcurves on both cylinders as certified B-splines, `Intersection` edges with cached 3D curves; tangency decided exactly; corpus cases with OCCT agreement. | 9 |
| M4 | Certificates | An interval bound of each branch's deviation from both surfaces, from the closed forms and core's outward-rounded `Interval`, stored as a `Certified` edge bound. | 3 |
| M5 | Controlled operators | The classical Euler set under ADR-0027's contracts: MVFS/KVFS, MEV/KEV, MEF/KEF, KEMR/MEKR, KFMRH/MFKRH, SEMV/JEKV, with radial rings, pcurve restriction and certificates; `MergeVertices` (ADR-0001); transactions as discarded candidate copies; local checks in checked builds; inverse-pair property tests; history events. ADR-0005's splitting edits join its benchmark. | 10 |
| M6 | Face clipping and imprinting | Branch against loop in parameter space (line and conic against line and circle pcurves, exactly decided); imprinting with M5. | 6 |
| M7 | Classification | Fragment connectivity; seeding by up to 16 deterministic rays per region, decided exactly for plane and cylinder faces; then the bounded generalised winding number (exact solid angles for planar faces, tessellated cylinders with a deviation bound), accepted only beyond its bound; then `Unresolved`; propagation; the cross-check table (ADR-0013's acceptance). | 10 |
| M8 | Assembly, merging, History | Keep rules; the ADR-0001 merge policy for vertices and edges; regularisation; a History with input ordinals, split and merge groups and Generated provenance (extending `HistoryEntry`). | 6 |
| M9 | Tolerance variants | Loosened inputs from 1e-4 to 1e-2 mm; growth logged; the merge-chain test. | 3 |
| M10 | Parallelism and determinism | Pair finding and surface intersection run in parallel through `ParallelFor`; the thread-count identity check; the 5× speed-up on a many-face fixture (a box with a grid of cylindrical holes); the G2 comparison job. | 4 |
| M11 | OCCT differential tests | Volumes, areas and face counts against `BRepAlgoAPI` on the same fixtures. | 2 |

About 65 working days in sequence. M0 comes first, and M1 needs none of it. M1, M2 and M5 have no dependencies on each other and can go in parallel; M3 needs M2, M3b needs M3 and M4, M6 needs M3 and M5, M7 needs M6, M8 needs M7, and M9 to M11 follow M8.

## 4. Decisions (2026-10-20, Onur Tuncer)

| # | Question | Decision |
| --- | --- | --- |
| D1 | Snapping (RS-0001 G3, G4) or nominal geometry (O4)? | **Nominal first.** Decide exactly on the input doubles, growing and logging a tolerance where a result needs it. Snapping comes as a `Context` option after RS-0001's go/no-go |
| D2 | The pair tolerance: the larger of the two (ADR-0001) or the sum (RS-0001 O1)? | **The larger**, as ADR-0001 states; RS-0001's O1 is answered |
| D3 | How are the prototype's edges certified? | **Interval bounds for the analytic branches now** (M4), stored as `Certified` |
| D4 | What are "loosened inputs"? | **Two fixture families**: loosened edges and vertices (raised tolerances, geometry perturbed within them) for ADR-0001, and surfaces moved within RS-0001's admissible placements for the corpus |
| D5 | The 1e-2 mm variant equals the ceiling | **That variant raises the ceiling in its `Context`**, as ADR-0001 says a caller does |
| D6 | An unresolved classification | **Fails the operation with a distinct diagnostic**, naming the fragment and the stage, and never guesses; a test asserts it does not happen on the Tier A corpus |
| D7 | ADRs 0010, 0013, 0027, 0029 and 0030 | **Reviewed and accepted before the work they govern starts** (M0), each narrowed to what the prototype builds |
## 5. Risks

- **Face clipping** in parameter space (M6) is where most Boolean bugs live: tangential contacts at loop boundaries, branches through vertices. Exact predicates help; the corpus needs cases for each.
- **History with two inputs** (M8) changes `HistoryEntry`. That is a public type, and `Commit` builds it today from one source.
- **The 5× speed-up** needs enough independent work: a box with a grid of 100 to 400 cylindrical holes is the planned fixture.
- **Full validation at every `Commit`** (ADR-0005) makes large results slow to publish. That is acceptable for the prototype; incremental validation is a Phase 1 item.
