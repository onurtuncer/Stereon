# ADR-0044: B-rep to tri-dexel conversion with bounded endpoints and explicit unresolved rays

- **Status:** Proposed
- **Date:** 2026-10-28
- **Owner:** Onur Tuncer
- **Phase:** 3–5 (staged; see Verification)
- **Related:** ADR-0001 (tolerances), ADR-0002 (robustness), ADR-0003 (capability tiers), ADR-0006 (errors), ADR-0007 (threading and determinism), ADR-0009 (pcurves), ADR-0010 (SSI), ADR-0013 (classification), ADR-0018 (validity), ADR-0021 (shape derivatives), ADR-0028 (purpose-specific tolerances), ADR-0029 (intersection outcomes), ADR-0033 (modular SDK), ADR-0043 (benchmarks)

## Context

Machining simulation represents in-process stock as a tri-dexel field: three
families of parallel rays, along X, Y and Z of a grid frame, each ray storing
the sorted intervals where it lies inside material. Stock removal then becomes
interval subtraction per ray, and cutter engagement a query on the intervals
near the tool. The field is initialised from the part or stock B-rep, so the
kernel must convert a closed solid into tri-dexels. The conversion has to be
fast, because it runs on large grids (10⁶–10⁸ rays). It also has to be accurate
in a sense the caller can check, because errors at this stage turn into gouge
or leftover stock in the simulation.

Three accuracies are easily confused:

- **Endpoint error:** the uncertainty of an interval endpoint along its ray.
- **Grid resolution:** the spacing between rays. Computing each endpoint to
  1 µm does not make a field with 100 µm spacing accurate to 1 µm across rays.
- **Feature preservation:** whether material that lies between rays of all
  three families (a thin rib, a small boss, a narrow slot) is seen at all.
  Tri-dexels catch features that are thin in one direction through the other
  two families; features small in all three directions can still fall
  between rays.

What exists today: the Boolean classifier casts rays against plane and
cylinder faces (`Detail::RayMeets`, `Detail::FaceCrossingsOf` in
`src/boolean/Classify.cppm`), with signed crossings and a refusal when a
meeting is ambiguous or a graze. There is no ray intersection for cones,
spheres, tori or B-spline surfaces, and no batched or coherent ray path.
`SampleCurve` and `SampleSurface` check their chord at seven interior points
per segment, a display bound that is not certified (`src/mesh/Sampling.cppm`),
so they cannot certify endpoints.

## Options

1. **Tessellate, then cast rays against the triangles.** Fast and mature, but
   the endpoint error is the chordal error plus the triangle-hit error, the
   mesh must be watertight at shared edges or rays leak, and the endpoints
   carry no face parameters or exact normals. A certified chordal bound would
   be needed for every face, so the cost moves into the tessellator.
2. **Cast every ray directly against the B-rep.** Endpoints come from the
   true geometry and can be certified, but without culling and coherence each
   ray tests many faces. Too slow at 10⁷ rays.
3. **Hybrid: conservative bounds for culling, exact geometry for endpoints,
   rays in coherent tiles.** Bounding hierarchies over faces and Bézier
   patches discard candidates per tile of rays. Analytic faces are solved in
   closed form, B-spline faces by subdivision with certified exclusion and
   interval Newton, and every hit is trimmed against the face's loops.
4. **Implicit or signed-distance voxelisation.** Simple and parallel, but its
   intervals are only as accurate as its sampling of the distance field. It
   has the same unbounded-between-samples problem as option 1.

## Decision

Adopt option 3. B-rep to tri-dexel conversion is a kernel operation of its
own library, `Stereon.Dexel`, above topology (ADR-0033). Batched ray to
trimmed-face intersection becomes a first-class kernel capability, built
ahead of general surface–surface intersection, which it does not need.
Toolpath strategy, cutter models, swept volumes and stock update policy stay
above the kernel. The field type and its interval operations (subtract,
intersect, unite per ray) belong to this library, so that applications share
one validated representation.

### 1. Input and grid

- The input is a body that passes the `Solid` validity level (ADR-0018):
  closed, oriented and manifold. Other bodies fail with `InvalidInput` and
  the checker's findings.
- The grid is explicit: an orthonormal frame with its origin, the spacing
  along each axis and the extent in cells. Rays of family *a* pass through
  cell centres of the plane normal to axis *a*. There is no implicit default
  grid. Spacings are lengths in the model unit (ADR-0001) and are recorded
  in the result.
- A request may cover a sub-range of tiles. Machining updates then rebuild
  only the tiles that a tool's swept bounds touch.

### 2. Output

- Per ray: sorted, pairwise disjoint material intervals `[enter, exit]`.
  Each endpoint carries its parameter along the ray and a **certified error
  bound** in model units, which is a numerical error bound in ADR-0028's
  sense, not a tolerance.
- Per endpoint, at the level of detail the caller asks for: the face, the
  surface parameters `(u, v)`, the outward unit normal and the contact kind
  (`Crossing`, `OnEdge`, `OnVertex`; ADR-0029). A field without metadata
  stores endpoints only. Metadata goes in parallel arrays, so the field
  stays compact.
- Storage is per family in compressed rows: one offset per ray, endpoints
  as `double`, metadata indices beside them. The layout is part of the
  public field type and is documented with it.
- The result records the grid, the input body's identity, the effective
  numerical policy (ADR-0028), the unresolved rays and the sub-resolution
  findings (§4).

### 3. Endpoint geometry

- **Analytic faces** in closed form: plane, cylinder, cone and sphere by
  their linear or quadratic equation, and the torus by a quartic with
  certified root isolation (interval Sturm sequence or Bernstein subdivision),
  never by an unchecked closed-form quartic.
- **B-spline faces** by Bézier extraction, done once per conversion and
  cached in the conversion's workspace, never in the `Shape` (ADR-0005). A
  patch is discarded only when its bounds *prove* that it misses the ray.
  Otherwise it is subdivided until interval Newton proves a unique root in
  the box, which then gives the certified bound. A box that is not resolved
  within the work budget makes the ray unresolved. Like ADR-0010's
  completeness rule, a patch is never skipped silently.
- **Trimming.** A hit counts only where its `(u, v)` lies inside the face's
  loops, decided with the robust predicates of ADR-0002 against the pcurves
  (ADR-0009). When the uncertainty box of the hit straddles a loop, the hit
  is resolved on the 3D edge (it is then `OnEdge`) or the ray is unresolved.
- **User geometry.** Tier-1 user surfaces (ADR-0003) have no certified ray
  intersection: a body that has one fails with `NotSupported` naming the face.

### 4. Intervals, events and feature preservation

- **Pairing by orientation, not by parity.** Each crossing is signed by the
  face's material side, as in `FaceCrossing::Entering` today. A ray's
  crossings, sorted, must alternate enter and leave, starting with enter and
  ending with leave. Hits on a shared edge or vertex are merged through
  topology into one crossing, or into none when both sides of it are the
  same (a graze along an edge of a convex or concave corner). A ray whose
  sequence does not alternate after merging is unresolved. Hits are never
  dropped or invented to make the sequence alternate.
- **Tangency.** A ray tangent to a face either touches it without crossing,
  so the face adds no endpoint, or lies in it (a ray in a planar face, a
  ray along a cylinder's ruling). Each is decided by the face's exact or
  certified predicate. When the predicate does not settle it, the ray is
  unresolved.
- **Unresolved rays** are listed per family, with the reason (ambiguous
  trim, exhausted budget, non-alternating sequence, undecided tangency) and
  the faces involved. By default an unresolved ray fails the conversion with
  `IntersectionFailed` and that list (ADR-0029). A caller may instead ask for
  a field that marks such rays. That field reports itself incomplete and is a
  distinct result, never a complete success.
- **Sub-resolution features.** After conversion, a conservative check
  reports every face whose bounds fall inside one cell column of all three
  families and that no ray hit. These are features the grid may have
  missed. The caller chooses to refine the grid locally, accept them as a
  finding, or fail. The check reports. It never moves material.

### 5. Performance

- Rays are processed in tiles (for example 16 × 16 rays of one family). A
  tile is culled against a bounding hierarchy over faces, then over cached
  Bézier patches, by the tile's slab, so most faces are discarded once per
  tile rather than once per ray.
- Analytic solves run across the rays of a tile in SIMD lanes
  (`Stereon::Simd`), with the same scalar operations per lane as the scalar
  path, so lane results equal scalar results bit for bit, as in `EvalGrid`.
- Tiles run on the `Context`'s pool through `Detail::ParallelFor`, in
  chunks of fixed size, each writing its own buffer. The buffers are joined
  in tile order. Results are bit-identical for every scheduler and thread
  count (ADR-0007 G1).
- Tessellation may be used for culling only if its bounds are conservative
  (a certified chordal bound inflates every triangle's box). Endpoints
  always come from the exact geometry.

### 6. Derivatives

Following ADR-0021, an endpoint on a transverse crossing of a fixed face
branch has the derivative `dt/dp = -(n · S_p) / (n · d)`, where `S_p` is the
surface's parameter derivative at the hit, `n` the normal and `d` the ray
direction. The conversion declares this product in its capability catalogue,
and an initial implementation may declare it `NotSupported`. A change in a
ray's interval count, a hit at an edge or vertex, a tangency, or `n · d`
near zero are events. They are reported as `Nondifferentiable` or
`NearEvent`, never differentiated through.

## Consequences

- Ray to trimmed-face intersection, with certified endpoints and trims, is
  built for every analytic surface and for B-splines before general SSI.
  Point-in-solid classification (B4.4), mass-property estimates and
  minimum-distance queries can reuse it.
- The Boolean classifier's plane and cylinder ray code is generalised and
  moved to the shared ray library, so the two never disagree about a hit.
- New work: torus root isolation, Bézier patch caching and certified ray
  to patch roots, the bounding hierarchies, the tile scheduler, the field
  type and its interval operations, and the sub-resolution check.
- Memory grows with the metadata level. A 1000³ grid with endpoints only and
  two intervals per ray needs about 3 × 10⁶ rays × 4 doubles ≈ 96 MB, before
  offsets. Metadata can triple that. The field documents its per-endpoint
  cost.
- The default refuses to return a field with unresolved rays. An application
  that wants best-effort output must ask for the incomplete result and handle
  it explicitly.

## Verification

Before acceptance, on the reference toolchain:

1. **Analytic corpus.** Boxes, cylinders, cones, spheres, ring and horn
   tori, and Booleans of them, with intervals known in closed form. Every
   ray's intervals match, and every endpoint lies within its certified bound
   of the exact value.
2. **B-spline corpus.** Random rational and non-rational patches of degree
   2–7, trimmed. The reference is a high-precision computation (interval
   subdivision in extended precision). There are zero missed or spurious
   crossings, and every endpoint is within its bound.
3. **Adversarial rays.** Rays through edges and vertices, along seams, at
   poles and apexes, tangent to cylinders, spheres and tori, lying in planar
   faces, and grazing concave and convex corners. Each ray is correct or
   unresolved, and none is silently wrong.
4. **Feature preservation.** Ribs, bosses and slots below the grid spacing
   in one, two and three directions. Features thin in one or two directions
   appear in another family, and features small in all three are reported
   by the sub-resolution check.
5. **Volume.** For each family, the interval lengths times the cell area
   converge to the analytic or `BRepGProp` volume as the spacing shrinks, at
   the rate the cell size predicts.
6. **Determinism.** Fields are bit-identical between the serial scheduler
   and parallel runs at 1, 2, 4 and 8 threads.
7. **Derivatives.** Endpoint derivatives with respect to primitive
   parameters (radius, height, pole position) agree with finite differences
   away from events, and events are reported where the interval count changes.
8. **Benchmark B5.7.** Same grids, same rays, same body, against OCCT's
   `IntCurvesFace_ShapeIntersector` on each ray, and for information against
   ray casting on a `BRepMesh` tessellation at matching chordal tolerance.
   Report rays per second at 1 and 8 threads, peak memory, endpoint agreement,
   unresolved and failed rays, and cold and warm runs, under ADR-0043's
   evidence rules. The target is ≥ 3× OCCT at 8 threads at equal endpoint
   accuracy, with zero silently wrong rays.

Stage the work: analytic faces and plane–cylinder bodies first (the
Boolean prototype's domain), then cones, spheres and tori, then trimmed
B-spline faces. Each stage passes items 1–6 for its domain before the next
one starts.
