# ADR-0048: A native ellipsoid surface

- **Status:** Proposed
- **Date:** 2026-11-05
- **Owner:** Onur Tuncer
- **Phase:** 2
- **Related:** [ADR-0001](0001-tolerance-model.md), [ADR-0003](0003-geometry-polymorphism.md), [ADR-0009](0009-parameter-space-representation.md), [ADR-0010](0010-ssi-strategy.md), [ADR-0021](0021-shape-derivatives.md), [ADR-0029](0029-intersection-and-classification-outcomes.md), [ADR-0047](0047-numerical-operation-completion-and-publication.md)

## Context

Stereon's analytic surfaces are the plane, cylinder, cone, sphere and torus. OCCT has the same five, and no ellipsoid either. In OCCT an ellipsoid is made by scaling a sphere unevenly with `BRepBuilderAPI_GTransform`, which turns its faces into B-spline surfaces. Every intersection and Boolean involving it then goes through the generic B-spline algorithms (`IntPatch`, `IntWalk`), with their marching tolerances.

An ellipsoid is a quadric:

> Σᵢ ((p − c) · eᵢ)² / aᵢ² = 1, with the semi-axes aᵢ along the orthonormal frame eᵢ.

Multiplying through by (a₁a₂a₃)² gives a polynomial in the input doubles. So everything Stereon decides exactly for the sphere can be decided exactly for the ellipsoid: a line's count of crossings, its touch, and a curve lying in the surface. The machinery is already there:
- `SignOfFormula` decides the sign of such a polynomial;
- the Bernstein subdivision of a curve's composed implicit equation (`Bernstein.hpp`) finds and proves its roots;
- the parametric–implicit tracer (`ImplicitTrace.cpp`) handles a quadric implicit side.

A native ellipsoid is therefore a capability OCCT does not have, at a cost bounded by machinery that exists.

What does not carry over from the sphere is the closest-point projection. A sphere's nearest point lies on the ray from its centre, but an ellipsoid's comes from a scalar root-finding problem that has no closed form.

## Options

1. **Represent ellipsoids as rational B-spline surfaces only**, as OCCT does. No new type, but every operation waits for the generic B-spline path (ADR-0010's Phase 2), and no decision about an ellipsoid can be exact.
2. **A general quadric type** (any second-degree surface by its matrix). One type would cover ellipsoids, hyperboloids and paraboloids, but it has no single natural parametrisation, chart or pole structure to give pcurves and topology (ADR-0009) a stable meaning.
3. **A native `Ellipsoid` with the sphere's parametrisation, chart and poles**, three semi-axes in place of one radius. Other quadrics could follow as their own types if needed.

## Decision

Adopt option 3.

### 1. The type

```
Center + cos v (A cos u XAxis + B sin u YAxis) + C sin v Axis(e)
```

- u ∈ [0, 2π] is the longitude and v ∈ [−π/2, π/2] the latitude, as for `BasicSphere`.
- `XAxis` and `YAxis` are orthonormal. `Axis` is their cross product when `RightHanded`, and its negative otherwise.
- A, B and C are the semi-axes along `XAxis`, `YAxis` and `Axis`.

`MakeEllipsoid` returns InvalidInput unless every value is finite, every semi-axis exceeds the Context's resolution, and the frame is orthonormal within FrameAllowance.

Equal semi-axes are accepted and kept as an ellipsoid: nothing is silently converted to a sphere. A handler may take a sphere's path when it proves A = B = C exactly.

Like the sphere's, the parametrisation is singular at the poles v = ±π/2. Their edges are ADR-0009's Degenerate edges, and the normal there is the limit ±`Axis`.

`Transformed` accepts a similarity, as for the other analytic types, and keeps the parameters. Uneven scaling of a sphere or an ellipsoid into an ellipsoid is a separate constructor, not `Transformed`. Its result has new axes (from the image of the frame and, in general, a symmetric eigenproblem) and a different parametrisation, so it must say so.

### 2. Projection

`Project(ellipsoid, p)` returns the parameters of a closest point. In the ellipsoid's frame, with p's coordinates yᵢ, the closest point is xᵢ = aᵢ² yᵢ / (aᵢ² + t), where t is the largest root of

> F(t) = Σᵢ (aᵢ yᵢ / (aᵢ² + t))² − 1.

F is strictly decreasing on (−min aᵢ², ∞), so t is found by bisection to the resolution of doubles: a fixed number of steps, deterministic, with no starting guess. That is Eberly's method ("Distance from a point to an ellipse, an ellipsoid, or a hyperellipsoid", 2013), whose 2D form `Ellipse.cppm` already cites. Where some yᵢ is zero the problem drops a dimension, as Eberly describes.

The closest point is not unique for points on the ellipsoid's medial set: the segment and ellipse inside it where two or more nearest points meet. `Project` then returns the one its bisection reaches, deterministically. Its documentation says that it promises a closest point, not continuity across that set (ADR-0047 §1).

Every caller that needs a point on the surface evaluates the returned parameters. None relies on the projection's own residual.

### 3. Intersections, in order

1. **Line × Ellipsoid.** In the frame, the line's quadratic has coefficients that are polynomials in the input doubles. The count, a touch and a line lying in the surface (impossible here) are decided exactly. The roots are placed from coefficients accurate to a few roundings (`QuadraticOf`), as for the sphere.
2. **B-spline, Circle, Ellipse × Ellipsoid.** The implicit equation composed with the curve's homogeneous poles is `EllipsoidForm`, the sphere's form with weighted squares. It is proved by the same subdivision.
3. **Plane × Ellipsoid.** A closed form: the section is an ellipse, a point or empty, decided by the sign of a polynomial in the inputs. The ellipse's centre and axes come from a 2 × 2 symmetric eigenproblem in the plane. The pcurve in the ellipsoid's chart is fitted and certified, since a planar section is not a conic in (u, v).
4. **Ellipsoid × the other analytic surfaces, and × Ellipsoid.** These go through the parametric–implicit trace, with `QuadricKind::Ellipsoid` as the implicit side. Closed forms for coaxial spheroids (sections that are circles) may follow, each behind its own exact decisions.
5. **Plane curves.** No new planar pair is needed: pcurves on an ellipsoid face are B-splines.

Each step lands with its own tests and lists what fails closed. Where a decision cannot be made, the result is `IntersectionFailed`, never an unproved empty result (ADR-0029).

### 4. Around the type

- **A solid primitive, `MakeEllipsoidSolid`:** one face, two Degenerate pole edges and a seam, mirroring the sphere's construction. Booleans then exercise the ellipsoid end to end. Winding numbers need its face integrals.
- **Derivatives (ADR-0021):** the parameters are the centre, the three semi-axes and the frame. Points of intersections on an ellipsoid get `F_z` and `F_p` like the sphere's.
- **Exchange:**
  - the Stereon archive stores the ellipsoid natively;
  - export to formats without an ellipsoid writes the exact rational B-spline, a sphere's rational patch scaled along its axes, whose image is exact;
  - import from them stays B-spline, as recognising an ellipsoid from a B-spline is not attempted.
- **Meshing, sampling, Python bindings, the viewer:** as for the sphere.

## Consequences

- An ellipsoid's intersections with lines, conics and B-spline curves, and with planes, are decided exactly or fail closed. They don't wait for the generic B-spline path, and OCCT reaches them only through that path.
- `BuiltinSurfaceTypes` grows by one. Every dispatch table, visitor, serialisation switch and binding that enumerates surfaces must handle it; the compiler finds the switches, and the dispatch tests find the tables.
- Projection is the one new numerical algorithm. Its non-uniqueness on the medial set is documented, not hidden.
- General affine images of ellipsoids, and other quadrics (hyperboloids, paraboloids), are out of scope.

## Verification

- **Exactness at the boundary:** decisions at 2^-40 gaps either side of tangency, for each handler.
- **Equal semi-axes as an oracle:** an ellipsoid with A = B = C must give the sphere's results, the same points to rounding with the same counts and contacts, for every curve type and surface pair.
- **Scaling an axis:** scaling a configuration by 2^k along one of the ellipsoid's axes, with axis-aligned inputs, maps a sphere's results onto the ellipsoid's exactly.
- **Adversarial grid:** the ellipsoid joins the grid of `tests/intersect/AdversarialTests.cpp`: validation as geometry, translation, scaling, reversal, far origins, seams and singleton domains.
- **Projection:**
  - closest-point optimality by sampling;
  - agreement with the sphere's projection when the semi-axes are equal;
  - the degenerate branches (some yᵢ zero, p on an axis, p at the centre);
  - determinism.
- **Booleans:** ellipsoid solids against boxes and cylinders through the integration checks of `tests/boolean/IntegrationTests.cpp`: closed shells, no duplicate edges, the volume identities, and the exact volume 4πABC/3.
- **OCCT differential:** the same configurations with OCCT's B-spline ellipsoid (a sphere through `BRepBuilderAPI_GTransform`), comparing points within OCCT's tolerance and recording where OCCT succeeds and Stereon fails closed, and the reverse.

Tier C owner review is required before implementation lands on `main`.
