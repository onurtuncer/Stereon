# Surface–surface relation cells: MR-01 to MR-04 for every ordered pair of surfaces

- **Date:** 2026-11-22
- **Branch:** `verify/relations-surface-surface`, on `main` 898ed69
- **Policy:** [ADR-0053](../adr/0053-metamorphic-testing.md), [inventory](intersection-metamorphic.json)
- **Toolchain:** GCC 16.2. The local disk was full during this work, so the evidence is the native CI of the pull request (Linux, Windows and macOS, GCC 16), not a local build.
- **Reviewer:** the v0.2 verification agent (Agent E). "Verified" means this documented review, pending the owner's review, as for the [planar](intersection-metamorphic-relations-planar.md), [curve–curve](intersection-metamorphic-relations-curve-curve.md) and [curve–surface](intersection-metamorphic-relations-curve-surface.md) cells.

The surfaces are `Plane`, `Cylinder`, `Cone`, `Sphere`, `Torus`, `Ellipsoid` and `BSplineSurface`: 28 unordered pairs, 49 ordered pairs, 196 cells. Before this branch 183 were `untriaged`. The other 13 were `MR implemented` by single-defect regressions: PR-244, PR-247, PR-230 and PR-230-pole, on Cone:Cone, Cone:Sphere, Sphere:Cone, Plane:Ellipsoid and Ellipsoid:Plane. Those regressions still run, and their defect records still cite them. The cells now cite the uniform tests below, which cover every pair alike.

All four relations apply to every pair. The service `Intersect(a, boxA, b, boxB, tolerance, ctx)` takes two surfaces in either order, so MR-02 is an exposed operation. Unlike curve–surface, no MR-02 cell here is an exception.

#### Domain

One transverse configuration per unordered pair, its answer known in closed form, never taken from the kernel:

- **The 15 pairs with a cone, a torus or an ellipsoid:** a case of [`corpus/curved-ssi/`](../../corpus/curved-ssi/README.md) at its `origin` placement and scale 1, each answered on `main` with no `may_fail_closed`, chosen to have at least one curve:

  | Family | Configuration | Curves |
  | --- | --- | --- |
  | cone-cone | coaxial-unequal | 2 |
  | cone-ellipsoid | both-nappes | 2 |
  | cone-torus | coaxial-two | 2 |
  | cylinder-cone | cone-axis-inside | 2 |
  | cylinder-ellipsoid | between-sides | 2 |
  | cylinder-torus | over-the-tube | 2 |
  | ellipsoid-ellipsoid | shifted-x-2 | 1 |
  | plane-cone | tilted-20 | 1 |
  | plane-ellipsoid | oblique-1p0 | 1 |
  | plane-torus | two-ovals | 2 |
  | sphere-cone | off-axis-loop | 1 |
  | sphere-ellipsoid | around-long-axis | 2 |
  | sphere-torus | across-the-tube | 2 |
  | torus-ellipsoid | through-the-tube | 2 |
  | torus-torus | shifted-1 | 2 |

  Each case's topology and witnesses come from `tools/generate_surface_corpus.py`'s closed-form arguments, as the corpus README states.
- **Two B-spline surfaces:** [`corpus/nurbs-ssi/`](../../corpus/nurbs-ssi/) `cylinder-cylinder/two-loops`, two rational quadratic cylinders of radii 1 and 2 with perpendicular axes, two loops; witnesses from `tools/generate_nurbs_corpus.py`.
- **The six pairs of planes, spheres and cylinders,** written in the test with their closed forms:
  - plane–plane: z = 0 over [-4, 4]^2 and x = 1/2 over y in [-3, 3]; the segment x = 1/2, z = 0.
  - plane–sphere: z = 1/2 and the unit sphere; the circle of radius sqrt(3/4).
  - plane–cylinder: the plane through (0, 0, 1/2) with axes x and (0, 3/5, 4/5) and the unit cylinder about z; the ellipse z = 1/2 + 4y/3.
  - sphere–sphere: unit spheres a unit apart along x; the circle x = 1/2 of radius sqrt(3/4).
  - sphere–cylinder: the unit sphere about (3/2, 0, 0) and the unit cylinder about z; z^2 = 3x - 9/4 on the cylinder, one loop over x in [3/4, 1], away from the sphere's poles.
  - cylinder–cylinder: the unit cylinder about z and the cylinder of radius 1/2 about x; two loops, x = ±sqrt(1 - y^2) with y^2 + z^2 = 1/4.
- **A B-spline surface against each analytic surface:** the bilinear B-spline that is exactly the plane z = 1/2 over [-6, 6]^2, against the plane x = 1/2, the cylinder of radius 2 about z, the cone of radius 1 and semi-angle 1/2, the sphere of radius 2, the torus (3, 1) and the ellipsoid (3, 2, 3/2). These are the plane sections of `SplineAnalyticTests.cpp`, now also under rotations.

**Claim.** The result has no overlap; it has the configuration's count of isolated points and of connected curves (`ConnectivityOf`, ends joined within twice the pair tolerance). Every branch must:
- be crossing, bounded and certified, with the cause `CertificateWidening` and its bound within the resolution;
- have its pcurves on both surfaces (`PcurveDeviation`).

Every one of 64 samples per branch, and every point, must lie on each analytic operand by that surface's closed-form distance: |h| for a plane; ||p - c| - R| for a sphere; |rho - R| for a cylinder; ||rho - |R + h tan a|| cos a for a cone; the tube distance for a torus; the algebraic value times the largest semi-axis over two for an ellipsoid. The bound is four times the branch's bound plus the tolerance. A B-spline operand has no closed-form distance here; its pcurve check and the other operand's closed form stand for it. Every witness must lie within 10^-3 (times the scale) of the result.

A refusal (`IntersectionFailed`, `CertificationFailed`, `NotSupported`, `NotImplemented`) fails the test. Each configuration is answered at its origin placement, so a refusal under a relation is a relation failure, never a pass.

**Out of scope:** tangencies, isolated contact points, overlaps and coincidences, near-degenerate configurations (the corpus's `may_fail_closed` cases), and seams or poles chosen to coincide with a branch. These categories are exercised, in their supported placements, by the corpus runners (`SurfaceCorpusTests.cpp`, `CorpusTests.cpp`) and the handler tests. They are not claimed by these cells.

#### Relations as the tests run them

The four tests are in `tests/intersect/SurfaceRelationTests.cpp`. Every claim is run in both operand orders, so each ordered pair has each relation.

- **MR-02:** "MR-02: surface-surface answers hold with the operands exchanged". The configuration as given, in both orders.
- **MR-01:** "MR-01: surface-surface answers move with exact rigid motions". Three exact rotations: the two cyclic permutations of the axes and a quarter turn about z. The first and last are composed with the translation (8, -4, 2). Points and centres move, unit frames turn, and the witnesses move with them. The translation rounds a few non-dyadic origins by about 10^-15, which is far inside every claim's margin for these transverse configurations.
- **MR-03:** "MR-03: surface-surface answers scale with the geometry and its tolerances". Scales 1/4 and 4, exact. All lengths are scaled: radii, a plane's box and a cylinder's or a cone's heights. Angles are not. The pair tolerance, the Context's resolution and its tolerance ceiling are scaled too (ADR-0053: the policy is transformed separately from the geometry).
- **MR-04:** "MR-04: surface-surface answers keep their place in equivalent surface charts". Each operand in turn is put in an equivalent chart, in both orders:
  - an analytic surface's frame is turned half about its own axis, X and Y negated, which leaves the point set and its normal unchanged. The chart becomes u' = u - pi modulo a turn (on a plane, (u, v)' = -(u, v)). A requested box is mapped the same way and kept within the base turn [0, 2 pi]; a whole turn stays [0, 2 pi]. A box outside the base turn is excluded from these cells, because of the finding below;
  - a B-spline's u is reversed: knots u' = a + b - u, rows of poles and weights reversed.

  The pcurves must lie on the recharted surface, so the reported parameters must be the remapped ones.

#### Replay

```sh
cmake --preset gcc16-debug
cmake --build --preset gcc16-debug --target stereon_intersect_tests
build/gcc16-debug/tests/stereon_intersect_tests "MR-0?: surface-surface*"
```

Every fixture is a committed corpus file or a literal in the test; there is no random seed. The tests are tagged `[intersect][surface-surface][metamorphic]`, not hidden, so they run in every pull request's CTest on Linux, Windows and macOS.

#### Finding: plane-cylinder boxes outside the base turn

The first CI run of this branch mapped a whole-turn box [0, 2 pi] to [-pi, pi] under MR-04. In all three runs (macOS and Windows `gcc16-debug`, Windows `gcc16-windows`: [macOS job](https://github.com/onurtuncer/Stereon/actions/runs/45059668981/job/137322549137), [Windows job](https://github.com/onurtuncer/Stereon/actions/runs/45059669001/job/137322549521)), every pair answered correctly except one: plane–cylinder, in both orders, **succeeded with half the ellipse**. The witness (0, 1, 11/6) lay 1.94 from the result.

- **Reproduction:** the plane through (0, 0, 1/2) with axes x and (0, 3/5, 4/5), box [-4, 4]^2. Take the unit cylinder about z with its frame turned half, X = (-1, 0, 0), Y = (0, -1, 0), and box u in [-pi, pi], v in [-4, 4]. Pass tolerance 1e-5 and the default Context. The answer is an ellipse arc over u in [0, pi] only. The closed form is the whole ellipse z = 1/2 + 4y/3, because [-pi, pi] is a whole turn of the same cylinder.
- **Cause, from reading the code:** `EllipseSection` in `src/intersect/PlaneCylinder.cpp` intersects the requested u range with the base turn [0, 2 pi) (`Intersection(Interval{-u0, CirclePeriod - u0}, ...)`) and drops the rest of the request without failing. `CircleSection` does the same (`Intersection(Domain(circle), boxB.X)`) and presumably has the same fault for a plane across the axis. This was not run.
- **Why it is a defect:** the request is a whole turn of the cylinder, and the answer is a successful, incomplete result: a wrong answer, not a refusal. Other handlers either answer such boxes (every other pair in this test did) or document a refusal: the ellipsoid's `NotImplemented` for a box outside its base turn (`docs/intersection.rst`).
- **Not repaired here:** this pull request changes no library code. The inventory records it as `AUDIT-20261122-plane-cylinder-periodic-box`, `fixture captured`. Its repair belongs to the SSI owner, and the owner decides whether such boxes are supported or refused.
- **Effect on the cells:** Plane:Cylinder/MR-04 and Cylinder:Plane/MR-04 stay `MR implemented` until that record is closed. The MR-04 test now keeps boxes within the base turn, so it no longer exercises this fault.

#### Green

At this branch's head `bbec12a` (main 898ed69 with test-only changes), the four tests passed in every native CTest job of #329: [Linux `gcc16-debug`](https://github.com/onurtuncer/Stereon/actions/runs/45065020641/job/137339209208), [Linux `gcc16-asan`](https://github.com/onurtuncer/Stereon/actions/runs/45065020641/job/137339209139), [Windows `gcc16-debug`](https://github.com/onurtuncer/Stereon/actions/runs/45065020649/job/137339207223), [Windows `gcc16-windows`](https://github.com/onurtuncer/Stereon/actions/runs/45065020649/job/137339207144), [macOS `gcc16-debug`](https://github.com/onurtuncer/Stereon/actions/runs/45065020688/job/137339132173) and [coverage](https://github.com/onurtuncer/Stereon/actions/runs/45065020698/job/137339134895). A refusal fails each test, so every configuration answered under every relation, in both orders.

The first run, at `c766f89`, was red only on the plane–cylinder MR-04 runs described in the finding above. MR-01, MR-02 and MR-03 were green there too.

**Review.** 194 of the 196 cells are `verified` on this agent's review of the fixtures, oracles and these runs, within the domain stated above. Plane:Cylinder/MR-04 and Cylinder:Plane/MR-04 stay `MR implemented` until AUDIT-20261122-plane-cylinder-periodic-box is closed. No fault spot-check was run for these cells, because the disk was full and nothing could be built locally. The owner's review is pending.
