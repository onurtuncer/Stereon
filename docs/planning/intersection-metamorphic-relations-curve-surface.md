# Curve–surface relation cells: MR-01 to MR-04 for every curve against every surface

- **Date:** 2026-11-12
- **Branch:** `verify/relations-curve-surface`, stacked on the curve–curve relations branch (#275)
- **Policy:** [ADR-0053](../adr/0053-metamorphic-testing.md), [inventory](intersection-metamorphic.json)
- **Toolchain:** GCC 16.2 (MSYS2 UCRT64), preset `gcc16-debug`, Windows 11
- **Reviewer:** the v0.2 verification agent (Agent E). "Verified" and "exception reviewed" mean this documented review, pending the owner's review, as for the [planar](intersection-metamorphic-relations-planar.md) and [curve–curve](intersection-metamorphic-relations-curve-curve.md) cells.

There are 42 ordered pairs: `Line`, `Circle`, `Ellipse`, `BSplineCurve`, `Hyperbola` and `Parabola` against `Plane`, `Cylinder`, `Cone`, `Sphere`, `Torus`, `Ellipsoid` and `BSplineSurface`. That makes 168 cells:

| Cells | Count | State |
| --- | --- | --- |
| MR-01, MR-03, MR-04 on the 28 pairs with a handler (line, circle, ellipse, B-spline curve against every surface) | 84 | verified, new tests |
| MR-02 on all 42 pairs | 42 | exception reviewed: no exposed exchanged operation |
| MR-01, MR-03, MR-04 on hyperbola and parabola against every surface | 42 | left open at first (no handler); handlers since #292, qualified in the 2026-11-22 addendum below |

There was no curve–surface corpus to cite. The oracle is a closed-form fixture in `tests/intersect/CurveSurfaceRelationTests.cpp`.

#### Domain

- **Surfaces:** each surface in its canonical frame, with a point P at chart parameters away from seams and poles, and its unit normal n. P and n come from the surface's defining formula, evaluated in the test:
  - plane: P = (0.3, 0.2, 0);
  - cylinder: radius 1.5, u = 1, height 0.7;
  - cone: radius 1, semi-angle 0.5, u = 1, v = 1.5;
  - sphere: radius 2, u = 1, latitude 0.6;
  - torus: radii 3 and 1, u = 1, v = 0.7;
  - ellipsoid: semi-axes 3, 2 and 1.5, u = 1, latitude 0.5;
  - B-spline surface: a bilinear patch that is exactly the plane z = 0 over [-4, 4]^2.
- **Curves:** each curve kind passes through P transversally (tangent d = cos 0.4 n + sin 0.4 e, with e across n) at a known parameter:
  - a line from P - 2d: t* = 2;
  - a circle (radius 3) and an ellipse (3 by 2) through P at angle 0.6, their frames turned for it;
  - a cubic B-spline, polynomial and rational, that is exactly c(u) = P + L (u - 3/8) d + h (u - 3/8)^2 x: u* = 3/8. The rational one carries the linear weight 1 + u in numerator and weights, so it is the same curve.
  - a hyperbola (semi-axes 3 and 2) and a parabola (focal length 1), added in the 2026-11-22 addendum: axes x and d, through P at t* = 0.6 (the parabola's times the scale).

  The crossings are deliberately off the kernel's piece boundaries (quarter-turn arcs, the dyadic subdivision); see the observation below.
- **Claim:** a reported point within four resolutions of P with curve parameter t*, no overlap, and every reported point on the curve and on the surface. The surface distance is closed-form: |z| for the planes; |rho - R| for the cylinder; ||rho - |R + z tan a|| cos a for the cone, both nappes; ||p - c| - R| for the sphere; the tube distance for the torus; and the algebraic value times the least semi-axis over two for the ellipsoid, an estimate near it.
- **3D rational curves:** included. The rational cubic is exercised against every surface.
- **Out of scope:** tangencies and overlaps (exercised by the handler tests themselves), unbounded domains, periodic B-splines, and unbounded hyperbola and parabola domains (`NotSupported`). Hyperbolas and parabolas had no handler until #292.

#### Relations

- **MR-01:** "MR-01: curve-surface crossings move with exact rigid motions". The identity, the two cyclic axis permutations and a quarter turn about z, each with and without the translation (8, -4, 2).
- **MR-03:** "MR-03: curve-surface crossings scale with the geometry and its resolution". Scales 1/4, 1 and 4, with the resolution and tolerance ceiling scaled. A parabola's parameter is a length and scales: t* = 0.6 times the scale.
- **MR-04:** "MR-04: curve-surface crossings keep their place in an equivalent curve chart". Under each rotation, the curve in its original chart and in an equivalent one:
  - a line's origin at P with its direction doubled, t' = (t - 2) / 2;
  - a circle's frame turned a quarter, t' = t - pi/2;
  - an ellipse's turned half, t' = t - pi;
  - a hyperbola's or a parabola's y axis negated, t' = -t;
  - a B-spline's knots mapped by u' = 3 - 2u (reversed, on [1, 3]).

  The point and parameter must follow.
- **MR-02:** exception reviewed. ADR-0053 section 1: "MR-02 applies only if reversing the operands is an exposed and semantically compatible operation". The service is `Intersect(curve, domain, surface, box, ctx)` only; there is no surface-curve operation to exchange the operands into. The same holds for every one of the 42 pairs.

#### Replay

```sh
cmake --build --preset gcc16-debug --target stereon_intersect_tests
build/gcc16-debug/tests/stereon_intersect_tests "MR-01: curve-surface*"
build/gcc16-debug/tests/stereon_intersect_tests "MR-03: curve-surface*"
build/gcc16-debug/tests/stereon_intersect_tests "MR-04: curve-surface*"
```

The fixture is deterministic. The tests run in every pull request.

#### Green

On this branch (main 227d79a with test-only changes):

| Test | Assertions | Runs answered | Refused (fail-closed) | No handler (hyperbola, parabola) |
| --- | --- | --- | --- | --- |
| MR-01 | 1532 | 296 | 0 | 96 |
| MR-03 | 567 | 111 | 0 | 36 |
| MR-04 | 1532 | 296 | 0 | 96 |

All three pass.

#### Observation: roots on a piece boundary refuse

A first version of the fixture put the crossings at a circle's or ellipse's angle 0 and at a B-spline's u = 1/2. Those are the kernel's piece boundaries: a quarter-turn arc cut, counted from the domain's lower end, and the first Bezier split. Every such run was refused with `IntersectionFailed`, "unresolved root existence or multiplicity", for circles, ellipses and both B-spline curves against most surfaces. This is fail-closed, never a wrong answer. The ledger documents it for a B-spline root at a knot, but not for a conic root exactly at an arc cut. It is a completeness gap for curves crossing exactly there, not a relation failure, and it is reported to the curve owner.

#### Spot-checks

Each fault was applied alone to a working tree with no committed change to `src/`, rebuilt, and run twice; the source was then restored.

- **Segment origin dropped.** The fault: `ParameterAt` in `Bernstein.hpp` maps a segment's local parameter back without the segment's start. It is killed in both runs by all three tests:
  - in line × cone and line × torus, whose segments start at -10;
  - additionally, **only by MR-04**, in every spline and rational pair (12 pairs). Their original chart starts at 0, which masks the fault; the equivalent chart on [1, 3] exposes it.
- **Conic arc offset dropped.** The fault: the second and later quarter-turn arcs start from 0 instead of the domain's lower end. The test process crashed (exit 0xC0000409) in every run. By `tools/mutate.py`'s rules that is inconclusive, not a kill, so it is not claimed as evidence.

#### Audit hardening, 2026-11-12

Every supported fixture must now return a result. `IntersectionFailed`,
`CertificationFailed`, and `NotImplemented` are failures for a supported pair,
not a warning counter. Only hyperbolas/parabolas against analytic surfaces may return
`NotImplemented`; the generic B-spline-surface handler must answer those conics
too. The six hyperbola/parabola versus B-spline relation cells remain open
pending their own qualification; they must not be described as missing handlers. This prevents lost handlers or new refusals from
preserving a green relation test.

Native audit validation (GCC 16.2 debug): the six curve-curve and curve-surface
relation tests pass 14,388 assertions. A controlled replacement of every
curve-surface answer by `IntersectionFailed` fails all three curve-surface
relation tests (exit 42); the unmodified source is restored afterwards.

#### Addendum, 2026-11-16: hyperbola and parabola handlers (#292)

#292 adds handlers for hyperbolas and parabolas against the analytic
surfaces (`OpenConicSurface`). `CurveSurfaceRelationTests.cpp` no longer
expects `NotImplemented` for them: every pair must now answer. On #292's
branch (main b292fd3 merged; Windows, GCC 16.2, `gcc16-debug`):

| Test | Assertions | Runs answered |
| --- | --- | --- |
| MR-01 | 2461 | 392 |
| MR-03 | 919 | 147 |
| MR-04 | 2461 | 392 |

All three pass. The runs that previously had no handler (96, 36 and 96) all
answer. The 36 hyperbola/parabola cells against analytic surfaces have
therefore been exercised, but they have not been reclassified here. With the
six against the B-spline surface, they await the same review as the 84 cells
above before the table at the top changes.

## Hyperbola and parabola cells, 2026-11-22

- **Branch:** `verify/relations-open-conics`, on `main` 898ed69
- **Reviewer:** the v0.2 verification agent (Agent E), as above; the owner's review is pending.
- **Cells:** MR-01, MR-03 and MR-04 of `Hyperbola` and `Parabola` against all seven surfaces, 42 cells. Their MR-02 cells stay exception reviewed, as for every curve–surface pair.

#### Fixture change

Until this branch the hyperbola and the parabola passed through P at their vertex, t* = 0. MR-04's chart for both negates the y axis, t' = -t, which leaves t* = 0 where it is. So the relation's parameter map was never checked at the crossing: a handler that reported the crossing's parameter with the wrong sign would have passed. The curves now pass through P away from the vertex:

- the hyperbola, semi-axes 3 and 2, axes x (across d) and d, centre P - 3 cosh(0.6) x - 2 sinh(0.6) d: t* = 0.6, and -0.6 in the equivalent chart;
- the parabola, focal length f equal to the scale, axes x and d, vertex P - t*^2/(4f) x - t* d with t* = 0.6 times the scale: t* and -t* in the equivalent chart.

The tangent at P is 3 sinh(0.6) x + 2 cosh(0.6) d, or t*/(2f) x + d = 0.3 x + d. Its component along the surface normal is at least 2.37 cos 0.4 - 1.91 sin 0.4 > 1.4 for the hyperbola and cos 0.4 - 0.3 sin 0.4 > 0.8 for the parabola, so each crossing is transverse. t* = 0.6 is off the hyperbola's unit stretches (cut at -1, 0 and 1 over [-2, 2]) and off their dyadic subdivision, and off the parabola's single piece [-4, 4] and its subdivision. Every other claim and tolerance of the test is unchanged.

#### Evidence

- **Tests:** the three cited tests, which run in every pull request (`[intersect][curve-surface][metamorphic]`, not hidden). Since the audit hardening of 2026-11-12 and #292, every fixture must answer: a refusal or `NotImplemented` fails the test.
- **Oracle:** the closed form above; the kernel is used only for the answer under test.
- **Green:** see "Green, hyperbolas and parabolas" below.

#### Green, hyperbolas and parabolas

At this branch's head `8df63a0` (main 898ed69 with this test change), the three cited tests passed in every native CTest job of the pull request (#324): [Linux `gcc16-debug`](https://github.com/onurtuncer/Stereon/actions/runs/45064765243/job/137338360038), [Linux `gcc16-asan`](https://github.com/onurtuncer/Stereon/actions/runs/45064765243/job/137338360152), [Windows `gcc16-debug`](https://github.com/onurtuncer/Stereon/actions/runs/45064765179/job/137338337967), [Windows `gcc16-windows`](https://github.com/onurtuncer/Stereon/actions/runs/45064765179/job/137338338354), [macOS `gcc16-debug`](https://github.com/onurtuncer/Stereon/actions/runs/45064765367/job/137338338048) and [coverage](https://github.com/onurtuncer/Stereon/actions/runs/45064765262/job/137338362846). Every fixture must answer, so each green run means every hyperbola and parabola pair answered with the crossing at t* (or -t* in the equivalent chart) under every motion, scale and chart.

**Review.** The 42 cells are `verified` on this agent's review of the fixture, the oracle and these runs. The limits are those of the curve–surface domain above: transverse crossings away from piece boundaries, bounded domains. No fault spot-check was run for these cells, because the disk was full and nothing could be built locally. The owner's review is pending.
