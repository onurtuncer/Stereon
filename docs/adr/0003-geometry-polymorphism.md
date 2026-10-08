# ADR-0003: Geometry polymorphism via external polymorphism

- **Status:** Proposed
- **Date:** 2026-09-28
- **Owner:** Onur Tuncer
- **Phase:** 0

## Context

The kernel handles many curve and surface types (line, circle, ellipse, B-spline, plane, cylinder, cone, sphere, torus, NURBS, extrusion, revolution, offset) and must let users add their own, such as an analytic blade surface. Classic inheritance (OCCT's `Geom_Curve` hierarchy) couples every type to a base class, forces heap allocation and reference counting, and makes value semantics impossible. Hot loops (tessellation, surface–surface marching) evaluate geometry millions of times, so dispatch cost matters.

## Options

1. **Inheritance hierarchy with virtual functions.** Familiar, open to extension; heap-allocated, pointer-heavy, intrusive.
2. **Closed `std::variant` of all built-in types.** Fast, value-semantic; closed to user types, and every new type touches every visitor.
3. **External polymorphism (type erasure).** Geometry types are plain structs with free functions; a `Curve`/`Surface` wrapper erases them behind an internal concept/model pair with small-buffer storage.
4. **Hybrid of 2 and 3.** Type-erased wrappers for the open, public API; a closed variant fast path for the built-in types inside hot loops.

## Decision

Adopt option 4.

- Geometry types (`Line`, `Circle`, `BSplineCurve`, `Plane`, `Cylinder`, `BSplineSurface`, …) are aggregates with no base class. Behaviour is provided by free functions found by ADL: `Eval`, `Deriv`, `Domain`, `BoundingBox`, `Project`, `Reversed`, `Transformed`.
- C++ concepts `CurveGeometry`, `Curve2dGeometry` and `SurfaceGeometry` define the required operations. Any type satisfying the concept can be wrapped.
- `Curve`, `Curve2d` and `Surface` are copyable value types with 96-byte small-buffer storage; larger types (NURBS with many poles) store an immutable shared payload so copies stay cheap.
- Hot loops call `VisitBuiltin(g, f)`, which dispatches once over a closed `std::variant`-like table of built-in types (generated with expansion statements and pack indexing, ADR-0016) and falls back to virtual dispatch for user types.
- Binary operations (curve–curve, surface–surface intersection) dispatch through a table keyed on type-ID pairs, with analytic special cases registered for known pairs and a generic NURBS fallback.

## Consequences

- New geometry types plug in without changing the kernel; users can add analytic surfaces for their own domain.
- Value semantics and SBO remove most heap allocation for simple geometry.
- Two dispatch paths must be kept consistent; property tests compare them on every built-in type.
- Debuggers show erased types poorly; pretty-printers are required.

## Verification

Benchmark a tight evaluation loop (10⁸ calls) on line, circle and cubic B-spline. Pass: the erased `Curve` path is within 10% of direct calls when using `VisitBuiltin`, and within 2× for the pure virtual fallback.
