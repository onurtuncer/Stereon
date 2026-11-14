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
| MR-01, MR-03, MR-04 on hyperbola and parabola against every surface | 42 | left open: no handler (`NotImplemented`) |

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

  The crossings are deliberately off the kernel's piece boundaries (quarter-turn arcs, the dyadic subdivision); see the observation below.
- **Claim:** a reported point within four resolutions of P with curve parameter t*, no overlap, and every reported point on the curve and on the surface. The surface distance is closed-form: |z| for the planes; |rho - R| for the cylinder; ||rho - |R + z tan a|| cos a for the cone, both nappes; ||p - c| - R| for the sphere; the tube distance for the torus; and the algebraic value times the least semi-axis over two for the ellipsoid, an estimate near it.
- **3D rational curves:** included. The rational cubic is exercised against every surface.
- **Out of scope:** tangencies and overlaps (exercised by the handler tests themselves), unbounded domains, periodic B-splines, and hyperbolas and parabolas, which have no handler.

#### Relations

- **MR-01:** "MR-01: curve-surface crossings move with exact rigid motions". The identity, the two cyclic axis permutations and a quarter turn about z, each with and without the translation (8, -4, 2).
- **MR-03:** "MR-03: curve-surface crossings scale with the geometry and its resolution". Scales 1/4, 1 and 4, with the resolution and tolerance ceiling scaled. A parabola's parameter is a length, but its t* = 0 is fixed.
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
