# ADR-0003: Geometry polymorphism via external polymorphism

- **Status:** Proposed
- **Date:** 2026-09-28
- **Revised:** 2026-10-10 (capability tiers for user geometry; SBO size is a tuning constant)
- **Owner:** Onur Tuncer
- **Phase:** 0

## Context

The kernel handles many curve and surface types (line, circle, ellipse, B-spline, plane, cylinder, cone, sphere, torus, NURBS, extrusion, revolution, offset) and must let users add their own, such as an analytic blade surface. Classic inheritance (OCCT's `Geom_Curve` hierarchy) couples every type to a base class, forces heap allocation and reference counting, and makes value semantics impossible. Hot loops (tessellation, surface–surface marching) evaluate geometry millions of times, so dispatch cost matters.

Openness has a limit: a surface that can only be evaluated point by point cannot take part in a boolean that promises completeness and certified tolerances (ADR-0009, ADR-0010). The subdivision search needs conservative bounds over parameter boxes, and the certificates need interval evaluation. A user type must declare which of these it provides, and the kernel must refuse rather than guess.

## Options

1. **Inheritance hierarchy with virtual functions.** Familiar, open to extension; heap-allocated, pointer-heavy, intrusive.
2. **Closed `std::variant` of all built-in types.** Fast, value-semantic; closed to user types, and every new type touches every visitor.
3. **External polymorphism (type erasure).** Geometry types are plain structs with free functions; a `Curve`/`Surface` wrapper erases them behind an internal concept/model pair with small-buffer storage.
4. **Hybrid of 2 and 3.** Type-erased wrappers for the open, public API; a closed variant fast path for the built-in types inside hot loops.

## Decision

Adopt option 4.

- Geometry types (`Line`, `Circle`, `BSplineCurve`, `Plane`, `Cylinder`, `BSplineSurface`, …) are aggregates with no base class. Behaviour is provided by free functions found by ADL: `Eval`, `Deriv`, `Domain`, `BoundingBox`, `Project`, `Reversed`, `Transformed`.
- C++ concepts `CurveGeometry`, `Curve2dGeometry` and `SurfaceGeometry` define the required operations. Any type satisfying the concept can be wrapped.
- **Capability tiers.** What a wrapped type may be used for depends on which optional concepts it also satisfies. The wrapper records the tier at construction, and an operation that needs a higher tier returns `KernelError::NotSupported` naming the type and the missing capability. The kernel never assumes a lossless NURBS equivalent exists for a user type.

  | Tier | Concept | Operations required | Enables |
  | --- | --- | --- | --- |
  | 1 | `CurveGeometry` / `SurfaceGeometry` | `Eval`, `Deriv`, `Domain`, `BoundingBox`, `Project`, `Reversed`, `Transformed` | display tessellation, mass properties, local queries, `Spatial` edges (ADR-0009) |
  | 2 | `BoundedGeometry` | `Bound(g, paramBox)`: a conservative enclosure of the image of a parameter sub-box; `DerivBound(g, paramBox)`; for surfaces `NormalCone(g, paramBox)` | subdivision search and pruning (ADR-0010), robust proximity and classification (ADR-0013), certificates (ADR-0009) |
  | 3 | `IntersectableGeometry` | either certified analytic intersection handlers registered in the pair table, or `ToBSpline(g, tol) -> expected<BSplineSurface, KernelError>` with a certified deviation | participation as a face surface in booleans and blends |

  All built-in types are tier 3. Tier 2 is what makes a user surface usable in any search; tier 3 is what makes it usable in booleans.
- `Curve`, `Curve2d` and `Surface` are copyable value types with small-buffer storage; larger types (NURBS with many poles) store an immutable shared payload so copies stay cheap. The initial buffer size is 96 bytes; it is a tuning constant fixed by the Phase 0 benchmark, not part of the public API or ABI promise.
- Hot loops call `VisitBuiltin(g, f)`, which dispatches once over a closed `std::variant`-like table of built-in types (generated with expansion statements and pack indexing, ADR-0016) and falls back to virtual dispatch for user types.
- Binary operations (curve–curve, surface–surface intersection) dispatch through a table keyed on type-ID pairs, with analytic special cases registered for known pairs and a generic NURBS fallback.

## Consequences

- New geometry types plug in without changing the kernel; users can add analytic surfaces for their own domain, and the tier tells them exactly what they get for the functions they implement.
- Users who want their surface in booleans must supply bounds and either analytic intersections or a certified B-spline conversion; that is real work, but the alternative is silent wrong answers.
- Value semantics and SBO remove most heap allocation for simple geometry.
- Two dispatch paths must be kept consistent; property tests compare them on every built-in type.
- Debuggers show erased types poorly; pretty-printers are required.

## Verification

Benchmark a tight evaluation loop (10⁸ calls) on line, circle and cubic B-spline. Pass: the erased `Curve` path is within 10% of direct calls when using `VisitBuiltin`, and within 2× for the pure virtual fallback. The same benchmark measures copy cost for buffer sizes 64, 96 and 128 bytes against the actual sizes of the built-in types and fixes the constant.

A tier-1-only test surface wrapped and passed to a boolean must return `NotSupported`; the same surface with `BoundedGeometry` and `ToBSpline` added must take part in the ADR-0001 prototype booleans.
