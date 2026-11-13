# Curved Booleans: cones, spheres, tori and ellipsoids (design note)

**Owner:** Onur Tuncer
**Author:** Agent G (curved Booleans), for the v0.2 release team
**Opened:** 2026-11-13, on main at `6180e7e`
**Status:** design, for review. No code changes with this note. Revised 2026-11-13 with the owner's answers to section 6: every case below is required for v0.2, none is an exclusion.
**Ledger rows:** *Intersection to Boolean* (Remaining: cones, spheres, tori and ellipsoids in the Booleans; coincident curved faces) and *Native ellipsoid* (Remaining: Booleans on doubly curved faces).
**Owner decision it implements:** curved Booleans stay in v0.2 (2026-11-10).

## 1. Where the Boolean stands

The prototype's pipeline (`docs/planning/boolean-prototype-plan.md`, section 2) is the same for every surface: candidate face pairs by box (`Bodies.cpp`), the surface intersection over the faces' parameter boxes (`Stereon.Intersect`), each branch clipped to both faces (`Clip.cpp`), the networks imprinted (`Imprint.cpp`), regions seeded by rays and, failing those, by the bounded winding number (`Classify.cpp`, `Winding.cpp`), and the kept regions assembled (`Assemble.cpp`). The intersection side already answers every pair of planes, cylinders, cones, spheres, tori and ellipsoids (ledger *Surface–surface*). What refuses is the Boolean's own knowledge of surfaces, which stops at planes and cylinders:

| Site | What it refuses or assumes |
| --- | --- |
| `Winding.cpp` `ScaleOf` (246, 328) and `AccumulateFace` (769 region) | `NotSupported` for a face that is not a plane or a cylinder; the deviation bound is a cylinder's sagitta about one axis |
| `Classify.cpp` `MeetSurface` | a ray against anything but a plane or a cylinder is ambiguous, so every ray of a seed facing a sphere is ambiguous and the seed falls to the winding number, which refuses |
| `Classify.cpp` `Entering`, `CoincidentSide`, `PlaneNormal`, `OnFace` | the material side of a hit is the nearer root (true for a cylinder, not for a cone's two nappes or a torus); coincidence seeds only on planes |
| `Bodies.cpp` `BoxOf` | a face's box is the hull of its edges; a sphere's only edges are its meridian seam and two poles, so the hull is a flat box in the seam's plane and pairs are missed silently |
| `Bodies.cpp` 345, 407 (`BoundaryBranch`, `AddCoincidence`) | `NotSupported` for coincident curved faces; a boundary is mapped only into a plane |
| `Assemble.cpp` 311 (`Carry`) | a use carried onto another parameter must be a line pcurve |
| `Assemble.cpp` 490 (`CopyLoop`) | a use without a pcurve (unchanged by this design: every coedge has one, ADR-0050 decision 2) |
| `Winding.cpp` `ExtentOf` | the extent is the hull of the vertices; a torus solid has one vertex, so the first target would already be the resolution |
| `build/Primitives` | no cone (frustum) or torus solid to test with |

The clip and the imprint are already generic in the curve kinds that the intersection gives on these surfaces: lines of constant v (coaxial circles), circles and ellipses on planes, and the trace's certified cubic B-splines on both sides, cut at both charts' seams (`ImplicitTrace.cpp`). Seams are handled as on the cylinder (#214). What they do not handle is a branch that reaches a pole or an apex (section 3.4).

## 2. Principles carried over

- **Nominal geometry, fail closed** (ADR-0029, prototype decision D1). Every new decision is a certified sign or an `Unresolved` / `NotSupported` / `IntersectionFailed`; nothing is guessed, no surface moves.
- **Rays first, winding number second, `Unresolved` third** (ADR-0013, decisions on acceptance), with the fixed budget of sixteen deterministic directions.
- **Pair tests by the sum** (ADR-0050 decision 5): `Coincident(distance, tol(a), tol(b))` for every on-face and on-edge test the extension adds.
- **Exact edges where the closed form exists** (ADR-0009 section 8): coaxial circles and plane conics stay `Exact`; traced branches are `Trim` edges in each imprint and `Intersection` edges once assembled (`Assemble.cpp` `KindFor`, unchanged).

## 3. Module by module

### 3.1 Winding number (`Winding.cpp`)

The method stays: faces triangulated in their parameter spaces, every face taking the same points along a shared edge, cells split until each is separated from the point by more than the face can stray from its triangles, the unseparated cells counted into the bound. Only the bounds are surface-specific, and today they are a cylinder's. They become second-order bounds that hold for any C2 chart:

- **`SurfaceScale` gains the second derivatives.** For each surface, bounds valid over a cell: `|Su|`, `|Sv|` (first order, the Lipschitz constant used now), a least stretch, and `Huu >= |Suu|`, `Huv >= |Suv|`, `Hvv >= |Svv|`. Closed forms, evaluated on the cell's interval (`Interval` arithmetic on the cell's v, so a cell near a pole gets its own, smaller bounds in u):

  | Surface | `|Su|` | `|Sv|` | `Huu` | `Huv` | `Hvv` |
  | --- | --- | --- | --- | --- | --- |
  | cylinder (now) | r | 1 | r | 0 | 0 |
  | sphere | r cos v | r | r cos v | r sin v (abs) | r |
  | ellipsoid (a, b, c) | max(a, b) cos v | max(a, b, c) | max(a, b) cos v | max(a, b) | max(a, b, c) |
  | cone | rho(v) = R + v sin a | 1 | rho(v) | sin a | 0 |
  | torus (R, r) | R + r cos v | r | R + r cos v | r | r |

  The ellipsoid's column is a safe bound by its largest equatorial and overall semi-axis, not a tight one; tightening it is not needed for correctness.
- **The sagitta becomes the second-order remainder.** A chord or triangle of parameter extent (du, dv) strays from the surface by at most (Huu du^2 + 2 Huv du dv + Hvv dv^2) / 8 for a segment, and the same with the factor three the code uses now for a triangle. `Sagitta(radius, du)` and its three callers (`CrossingOf`, `DeviationOf`, `SplitAxisOf`) take the cell's `SurfaceScale` instead of a radius. For a cylinder the new bound equals today's (Huv = Hvv = 0), which is what makes slice 0 a refactor with no change of output (section 5).
- **The edge rate takes both parameter directions.** `PcurveRates` bounds |u'| only, because only u turns on a cylinder. It gains |v'|, and `RateOf` uses Huu u'^2 + 2 Huv |u' v'| + Hvv v'^2 + (|Su| + |Sv|) |p''|.
- **Poles.** At v = ±pi/2 on a sphere or an ellipsoid (and at a cone's apex) the least stretch in u is zero, which would make `CheckMeeting` accept any gap in u. The meeting check is split by direction: a gap in v is held to the least stretch in v (r, or c, or 1 on the cone), and a gap in u is accepted only on a Degenerate edge's line, where the whole segment is one point in space; there the polygon bridges the gap with corners on that line, which are one point in space, so their triangles have no area and the faces still close. A cell that touches a pole needs no special case beyond that: its u-stretch bound cos v goes to zero with the cell.
- **Two seams.** A torus face is periodic in both u and v. Its polygon is one loop of four seam uses (u-seam, v-seam, both reversed), already closed in the parameter plane; nothing changes but the scale.
- **Extent.** `ExtentOf` takes the hull of the faces' boxes (section 3.3), not only of the vertices.

This is the "tessellated with a known deviation bound" of ADR-0013's decisions on acceptance, with the bound now given for every elementary surface.

### 3.2 Ray classification (`Classify.cpp`)

A ray meets each surface where its implicit equation vanishes along the line. The pattern is `MeetCylinder`'s: coefficients enclosed in `Interval` arithmetic from the input doubles, a certainly negative discriminant a miss, a certainly positive one two crossings, anything else ambiguous (the next direction); the parameters then constructed in doubles and checked by `Refined`, Newton steps on a residual that is a distance to first order.

| Surface | Equation along o + s d | Residual for `Refined` | Material side at a hit |
| --- | --- | --- | --- |
| sphere | quadratic, \|w + s d\|^2 = r^2 | \|p - c\| - r | nearer root enters |
| ellipsoid | quadratic in the ellipsoid's frame, sum ((w + s d) . e_i / a_i)^2 = 1 | signed distance by Eberly's projection (`Ellipsoid.cpp`) | nearer root enters |
| cone | quadratic, both nappes; sin a and cos a enclosed as intervals of the double angle | distance in the meridian plane (`Cone.cppm`'s closest point) | sign of d . grad F at the hit, certified |
| torus | quartic | distance to the tube's centre circle less r | sign of d . grad F at the hit, certified |

- **Spheres and ellipsoids** are the cylinder's case over again: a closed convex quadric, two roots, the nearer one entering.
- **Cones** cannot use "nearer enters": a line through both nappes leaves one and enters the other. The material side of every hit is the sign of d . grad F at the refined hit, enclosed in intervals; an unsettled sign is an ambiguous ray. A hit on the other nappe is found and then dropped by `CrossesFace`, since its parameters lie outside the face's loops. A hit within tolerance of the apex is ambiguous (the apex is a vertex).
- **Tori** have a quartic. Rather than a new quartic solver, the ray uses the curve–surface service, `Intersect(Line, [0, ∞), Torus)` (`LineSurface.cpp`), whose roots are enclosed and whose tangent contacts fail closed; any failure, tangency or root it cannot separate is an ambiguous ray. Hits are then refined and their sides read from d . grad F as for the cone.
- **Coincident seeds.** `CoincidentSide` today tests planes only. It becomes: the seed lies on a face of the other body when its distance to that face's surface, by the surface's projection, is within the pair tolerance of the two faces (ADR-0050) and its projection is strictly inside the face's loops; OnSame or OnOpposite by the sign of the product of the two faces' outward normals at that point (Su x Sv, turned for a reversed face; at a pole the limit ±axis, as `Sphere.cppm` defines it). Only coincidences the imprint bounded reach this (section 3.6).
- **Point-in-face** (`FaceContains`, `Clip.cpp`) tests "on the boundary" with the edge nearest in parameters, measured in space. Near a pole the parameter metric is stretched in u, so the nearest in parameters can be farther in space than the true nearest. The test is kept as a filter that may only fail towards "on the boundary" (ambiguous), by measuring the distance in space to the edge's nearest point over a bracket of parameters around the parametric nearest; a point it misses is still caught by the parity rays, which decide exactly in the parameter plane.

### 3.3 Pairing and the faces' boxes (`Bodies.cpp`)

`BoxOf` must bound the face, not its edges: a face is the image of its trimmed parameter domain, which lies in its parameter box (`ParameterBoxOf`), so the face lies in the surface's image of that box. For each surface that image's box is enclosed by evaluating the surface over the parameter box in `Interval` arithmetic (with whole turns in u or v bounded by the surface's own `BoundingBox`), hulled with the edges' boxes as now. Without this a sphere and a box beside its far side are never paired and the Boolean returns a wrong solid with no error, which is the most dangerous failure in this note; slice 1 holds it with a fixture that would be missed.

### 3.4 Clipping and imprinting (`Clip.cpp`, `Imprint.cpp`)

- **Periodic charts and seams** already work: `NearestOn` places a periodic parameter in an edge's range, the imprint tells the two uses of a seam apart by where the piece's pcurve ends (`NearestOf`), and the intersection cuts every branch at both charts' seams. A torus has two seams; the same mechanism applies to v.
- **Poles and apexes.** A branch that reaches a pole of a sphere or an ellipsoid, or a cone's apex, arrives at the Degenerate edge at some u0. The face's loop must then pass through the pole vertex twice, with the Degenerate edge split at u0 into two Degenerate edges on the same vertex. No operator does this: `SplitEdgeMakeVertex` (SEMV) refuses Degenerate edges, rightly, since it makes a new vertex. ADR-0027 (decisions on acceptance, third bullet) admits non-Euler operators with the same contracts, so this is a new contracted `ShapeBuilder` operator, *split a Degenerate edge at a parameter, keeping its vertex*, with its inverse (join two Degenerate edges of one vertex on one face), not a new ADR. It is in `src/topo`, so it needs the coordinator's allocation. The intersection side also refuses such branches today: the trace ends with `IntersectionFailed` at a branch through a pole of the carrier's chart (`ImplicitTrace.cpp`), and closed-form circles on a sphere whose poles are off the common axis go to the trace. Until both land, a Boolean whose intersection passes through a pole or an apex fails closed with `IntersectionFailed` from the intersection, and the tests hold that; slice 6, required for v0.2, builds both.
- **Point contacts** (a sphere resting on a plane, two spheres touching): the intersection gives isolated points (`SphereTouch`, `AddAxisPoint`); nothing is imprinted, the seeds decide, and the regularised result keeps nothing of a contact of lower dimension, as for the existing edge contacts (#212).

### 3.5 Assembly (`Assemble.cpp`)

`Carry` puts a use of one body's edge on the parameter of the edge it is identified with. It accepts lines only. The generalisation keeps exactness where it exists and fails closed elsewhere:

- a `Line2d` (as now): any affine change of parameter;
- a `Circle2d` or `Ellipse2d`: a change s = ±t + c only (a shift of the angle and a reversal), as a turn of its frame by c and, reversed, its Y axis negated; the turned pcurve is checked against the original at the range's ends and middle within the edge's tolerance, since cos c and sin c are rounded;
- a non-rational or rational `BSplineCurve2d`: an affine map of its knots, the poles and weights unchanged (an exact reparametrisation up to the rounding of the knots, checked the same way);
- anything else: `NotSupported`, as now.

Both bodies' imprints of one piece of the intersection come from one branch range (`Bodies.cpp` `AddPiece`), so the curved case arises mostly for coincident boundaries (section 3.6).

### 3.6 Coincident curved faces (`Bodies.cpp`, `Classify.cpp`)

Coincident planes work by imprinting each face's boundary onto the other through `PCurveIn(plane, edge curve)`, then seeding OnSame/OnOpposite. For a curved face the boundary's pcurve on the other surface comes from the overlap's map instead: the intersection returns an `OverlapRegion` whose `AToB` is an affine map of the parameters for coaxial alike surfaces of revolution (`ElementaryPairs.cpp` `Coincident`) and for coaxial cylinders (`CylinderCylinder.cpp`). The edge's own pcurve, mapped by `AToB` (or its inverse), is the branch's pcurve on the other face: a line stays a line, a circle a circle (the map is an isometry of the parameter plane, up to reflection), a B-spline maps pole by pole. The mapped pcurve is then cut where it crosses the target chart's seam, and each piece shifted into the chart, as the intersection does for its own branches. The rest (clip, network, pairing, `CoincidentSide` of section 3.2) is the planar path.

Dependencies on the intersection side (Agent B's files), all fail-closed today:

- `ElementaryPairs.cpp` `Coincident` returns `NotImplemented` unless both boxes are whole charts. A face of a Boolean result is seldom a whole chart, so coincident spheres, cones and tori over part of their charts are needed (ledger *Surface–surface*, Remaining). The map does not change with the box; the boxes become the overlap's `OnA`, `OnB` and the Boolean clips by the loops, as for planes.
- Spheres whose poles differ have no affine map between their charts. The boundary's pcurve on the other sphere is then not the mapped pcurve but the edge's own curve placed on the other chart: an exact 3D circle or the edge's S(p(t)), its pcurve fitted and certified on the target sphere with the fitter and `CertifyDeviation` the trace already uses (cut at the target's seam, and at its poles by slice 6). The intersection side must also report the overlap of such spheres (today `NotImplemented`), with `SameOrientation` from the normals and no affine `AToB`; the Boolean then takes the fitted route.
- Two equal ellipsoids go to the trace, whose f is identically zero, so it fails. A closed form is needed: same centre, the same axes up to sign and order, equal semi-axes, each decided exactly; the map is the identity or a reflection of the chart.

### 3.7 Primitives (`src/build/Primitives.cpp`)

The tests need a cone frustum and a ring torus as solids, each built like the cylinder and the sphere through `ShapeBuilder` and checked at the Solid level: `MakeConeSolid` (two circles, a seam, a conical face and two disks; an apex form with a Degenerate apex edge later, for slice 6) and `MakeTorusSolid` (ring tori, `MajorRadius > MinorRadius`: one vertex, a u-seam and a v-seam, one face; horn tori, `MajorRadius == MinorRadius`, which the surface admits: the circle v = pi collapses to the centre, so the v-seam is placed there as a Degenerate edge at the centre vertex and the face is bounded like a sphere's, the centre handled as a pole by slice 6). These are construction additions with their own tests and, being outside `src/boolean`, need the coordinator's allocation.

## 4. Test oracles

Every slice is held to the existing invariants of `tests/boolean/IntegrationTests.cpp` (Solid-level validity but L05, no coincident edges or vertices, the volume identities through the test measure `tests/support/ShapeMeasures.hpp`) and, in addition, to oracles that do not use the kernel:

- **Membership by closed form.** A grid point's side in each operand is decided from the operand's defining inequalities (|p - c| < r, the box's slabs, the ellipsoid's quadratic form, a cone's and a torus's implicit inequalities within their caps), and its side in the result by the operation; the kernel's winding number of the result must agree wherever the point is farther than the operands' tolerances from every surface. This replaces the integration tests' use of `WindingNumberOf` on the operands, which is the kernel checking itself.
- **Closed-form volumes,** compared with the test measure within its own estimate and a relative 1e-9:

  | Configuration | Volume |
  | --- | --- |
  | sphere and plane (cap of height h) | pi h^2 (3r - h) / 3 |
  | two spheres (lens), centres d apart | pi (r1 + r2 - d)^2 (d^2 + 2 d (r1 + r2) - 3 (r1 - r2)^2) / (12 d) |
  | sphere less a coaxial cylinder of radius a (napkin ring) | pi H^3 / 6, H = 2 sqrt(r^2 - a^2) |
  | sphere and an octant box at its centre, poles off the planes | pi r^3 / 6 |
  | ellipsoid and a plane z = h across its c axis | pi a b (c - h)^2 (2c + h) / (3 c^2) |
  | cone frustum cut by a plane across its axis | pi H (R1^2 + R1 R2 + R2^2) / 3 |
  | ring torus and the half-space z < h, \|h\| < r | 2 pi R (r^2 acos(-h / r) + h sqrt(r^2 - h^2)) (Pappus) |
  | ring torus less a coaxial cylinder through its hole | Pappus on the remaining segment of the tube's disk |

- **Metamorphic relations** (ADR-0053): operand exchange for Union and Common (same volume and membership), and vol(A - B) - vol(B - A) = vol(A) - vol(B) for the Cuts; a rigid motion by an exact Pythagorean rotation and a dyadic translation (volume and membership carried); scaling by 2^k (volume by 8^k, exactly representable). The ellipsoid with equal semi-axes against the sphere (ADR-0048's own oracle): the same Boolean, the same volume.
- **Fail closed while unbuilt.** Every case of section 5 is required for v0.2 (owner, 2026-11-13). Until its slice lands, a configuration that needs it (a branch through a pole or an apex, a horn torus's centre, coincident spheres with different poles, partial-chart or ellipsoid coincidences) must fail with `NotSupported` or `IntersectionFailed`, and a test asserts the code; the slice that builds it turns that test into a closed-form one. In random properties such a refusal is `RC_PRE`, never a pass and never a failure.
- **Red before, green after.** Each slice's fixtures are first committed against the unmodified Boolean, where they fail with the `NotSupported` they name (or, for the box of section 3.3, with a wrong volume), and pass after.

## 5. Slices

Each slice is one draft PR, stacked, end to end through Union, Common and both Cuts.

| Slice | Scope | Files | Oracle |
| --- | --- | --- | --- |
| 0 | Second-order bounds in `Winding.cpp` for the plane and the cylinder; face boxes from the surface (`Bodies.cpp`); extent from face boxes | `Winding.cpp`, `Bodies.cpp` | every existing Boolean test unchanged, winding numbers of the existing fixtures bit-identical or within the old bound |
| 1 | **Spheres** against planes and spheres: ray quadratic, sphere bounds, poles in the meeting check | `Classify.cpp`, `Winding.cpp`, tests | cap, lens (centres on the pole axis: exact circles; off it: traced), octant, sphere inside a box, apart, touching at a point; a box beside the sphere's far side (the pairing fixture); a pole-crossing configuration failing closed |
| 2 | **Ellipsoids** against planes, spheres and ellipsoids | `Classify.cpp`, `Winding.cpp`, tests | ellipsoid cap; equal-semi-axes against slice 1 |
| 3 | Spheres and ellipsoids against **cylinders** | tests, any `Carry` case found | napkin ring; off-axis cylinders through the volume identities |
| 4 | **Cones** (frustums): ray with certified nappe sides, cone bounds, `MakeConeSolid` | `Classify.cpp`, `Winding.cpp`, `Primitives.cpp` | frustum by a plane; cone and sphere coaxial; a branch through the apex failing closed |
| 5 | **Tori** (ring and horn): ray through the curve–surface service, torus bounds, `MakeTorusSolid` | `Classify.cpp`, `Winding.cpp`, `Primitives.cpp` | ring and horn tori by a plane across the axis (Pappus); torus and coaxial cylinder; a horn torus cut through its centre failing closed until slice 6 |
| 6 | **Poles, apexes and the horn torus's centre in the branch**: the Degenerate-edge split operator (own PR), the imprint at a pole, the intersection's pole-crossing branches (closed-form circles through a sphere's or an ellipsoid's poles with pcurves split at the pole, and the trace's branches through a pole) | `src/topo`, `Imprint.cpp`, `Clip.cpp`, `ElementaryPairs.cpp`, `ImplicitTrace.cpp` | hemisphere by a plane through the pole axis, 2 pi r^3 / 3; an ellipsoid halved through its poles, 2 pi a b c / 3; a cone cut through its apex by a plane holding the axis; a horn torus halved through its centre |
| 7 | **Coincident curved faces**: cylinders, spheres (with the same or different poles), cones, tori, ellipsoids, over whole and partial charts | `Bodies.cpp`, `Classify.cpp`, `Assemble.cpp` (`Carry`), `ElementaryPairs.cpp` (partial charts, different poles, the ellipsoid's closed form) | A op A for each solid (A, A, empty); coaxial cylinders overlapping in height, pi r^2 (h1 + h2 - overlap); a hemisphere (from slice 1) against its sphere; a sphere against the same sphere with turned poles; an ellipsoid against a half of itself |

Slices 0 and 1 are the first implementation PR's scope; slice 0 may merge into it if its diff is small.

## 6. What accepted ADRs do not decide, and open questions

None of the slices needs a new ADR or a change to an accepted ADR's contract. These points are not decided by an accepted ADR; the owner's answers of 2026-11-13 are recorded with them:

1. **"Exact predicates" for curved ray hits.** ADR-0013's decisions on acceptance name a ray against a plane or a cylinder, decided with exact predicates; the prototype implements that as signs enclosed in `Interval` arithmetic from the input doubles, an unsettled sign making the ray ambiguous. This design reads it the same way for spheres, ellipsoids and cones (a cone's sin a and cos a enclosed from the double angle) and, for the torus's quartic, as the enclosed roots of the curve–surface service. The fail-closed path (next ray, winding number, `Unresolved`) is unchanged. **Decided by the owner, 2026-11-13:** for curved rays ADR-0013's "exact predicates" means interval-certified signs, a sign decided only when an enclosure proves it and the ray failing closed otherwise. No ADR change.
2. **Splitting a Degenerate edge** (section 3.4) is a new contracted operator under ADR-0027's third decision on acceptance, not a new ADR. *Allocated by the coordinator (2026-11-13) to this stream, as its own small PR.*
3. **ADR-0048 is Proposed.** The ellipsoid's Booleans (slice 2) build on a proposed ADR, as the ellipsoid's surface, intersections and primitive already do. *Owner, 2026-11-13: ADR-0048 is being accepted; the coordinator changes its status in a separate documentation PR.*
4. **The v0.2 supported domain for curved Booleans.** **Decided by the owner, 2026-11-13: no exclusions.** v0.2 must solve every case: branches through a pole or an apex, coincident curved faces over part of a chart, coincident spheres whose poles differ, coincident ellipsoids and horn tori, as slices 5 to 7 plan them. Each may fail closed while it is unbuilt, never silently wrong.
5. **Tolerance growth of declared tangencies.** A coaxial tangency declared within the pair tolerance gives a branch whose bound is that tolerance; the imprint makes its edge at that bound without a `ToleranceChange` (ledger *0.2-T2*, Remaining). Curved Booleans make this case common (a sphere inscribed in a cylinder, cones touching). This is the open question already recorded there, not a new one.
6. **Cross-stream files.** Slices 4 to 7 touch `src/build/Primitives.cpp` and `src/topo` (both allocated to this stream, 2026-11-13) and `src/intersect/ElementaryPairs.cpp` and `ImplicitTrace.cpp` (free, to be announced to the coordinator before work starts); slices 0 to 3 stay in `src/boolean` and `tests/boolean`.

## 7. References

- ADR-0009 (edge kinds, Degenerate edges), ADR-0013 (classification), ADR-0027 (operators), ADR-0029 (outcomes and side transitions), ADR-0048 (ellipsoid), ADR-0050 (tolerant topology), ADR-0053 (metamorphic testing), ADR-0054 (events, proposed).
- `docs/planning/boolean-prototype-plan.md`, `docs/planning/v0.2-acceptance-ledger.md`.
- J. Van Oosterom and J. Strackee, "The solid angle of a plane triangle", *IEEE Trans. Biomed. Eng.* 30 (1983); A. Jacobson, L. Kavan and O. Sorkine-Hornung, "Robust inside-outside segmentation using generalized winding numbers", *ACM TOG* 32 (2013); C. Hoffmann, *Geometric and Solid Modeling* (1989), section 7.2.
