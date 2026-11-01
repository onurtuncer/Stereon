# ADR-0003: Geometry polymorphism via external polymorphism

- **Status:** Accepted
- **Date:** 2026-09-28
- **Revised:** 2026-10-10 (capability tiers for user geometry; SBO size is a tuning constant)
- **Revised:** 2026-10-20 (the virtual fallback's target is an absolute overhead per call, not a ratio; one buffer size per wrapper kind; tiers 2 and 3 arrive with the intersection and boolean work)
- **Accepted:** 2026-10-20 by Onur Tuncer: `VisitBuiltin` is within 10% of direct calls, and the virtual fallback adds about 3.4 ns per call, within the revised target (see *Results*)
- **Revised:** 2026-11-01 by Onur Tuncer, after the review of that day: `VisitBuiltin`'s target is measured against the same loop on a local copy of the geometry, the ratio to the direct call reported beside it (*Verification*); the binary dispatch is described as built (compile-time tables, no fallback yet); the tier check in the Booleans is in place (`tests/boolean/TierTests.cpp`)
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
  | 3 | `IntersectableGeometry` | either certified analytic intersection handlers in the pair tables, or `ToBSpline(g, tol) -> expected<BSplineSurface, KernelError>` with a certified deviation | participation as a face surface in booleans and blends |

  All built-in types are tier 3. Tier 2 is what makes a user surface usable in any search; tier 3 is what makes it usable in booleans.

  Phase 0 builds tier 1. The `BoundedGeometry` and `IntersectableGeometry` concepts are defined with the subdivision search and the booleans that use them (ADR-0010, ADR-0013). Until then an operation that would need them returns `NotSupported` for a user type, as the intersection dispatch table already does.

  Design-parameter derivatives are a separate optional capability, orthogonal to the tiers (ADR-0021).
- `Curve`, `Curve2d` and `Surface` are copyable value types with small-buffer storage; larger types (NURBS with many poles) store an immutable shared payload so copies stay cheap. Each wrapper kind has its own buffer, sized to its largest built-in type plus the model's table pointer, so that every built-in analytic type is stored inline: 96 bytes for `Curve` (`Ellipse`, `Hyperbola`), 72 for `Curve2d` (`Ellipse2d`, `Hyperbola2d`) and 104 for `Surface` (`Cone`, `Torus`). The sizes are tuning constants, not part of the public API or ABI promise. The copy benchmark first set one 96-byte buffer, then 128 when `Cone` and `Torus` arrived (a shared copy costs about three times an inline one); topology then made the wrapper's size a memory cost per edge and per coedge (ADR-0004, B0.4), which fixed the per-kind sizes (`docs/geometry.rst`, *The buffer size, revisited*).
- Hot loops call `VisitBuiltin(g, f)`, which dispatches once over a closed `std::variant`-like table of built-in types (generated with expansion statements and pack indexing, ADR-0016) and falls back to virtual dispatch for user types.
- Binary operations (curve–curve, plane curve, curve–surface and surface–surface intersection) dispatch through tables indexed by the pair of built-in type indices, generated at compile time from the built-in type lists. They are constant data, not a registry, so there is no global mutable state (ADR-0007). A pair with a closed-form method has its analytic handler; every other built-in pair has a stub that returns `NotImplemented` naming both types, and a type that is not built in gets `NotSupported`. The generic NURBS path that will replace the stubs belongs to ADR-0010 and is not built yet (`src/intersect/Dispatch.cppm`).

## Consequences

- New geometry types plug in without changing the kernel; users can add analytic surfaces for their own domain, and the tier tells them exactly what they get for the functions they implement.
- Users who want their surface in booleans must supply bounds and either analytic intersections or a certified B-spline conversion; that is real work, but the alternative is silent wrong answers.
- Value semantics and SBO remove most heap allocation for simple geometry.
- Two dispatch paths must be kept consistent; property tests compare them on every built-in type.
- Debuggers show erased types poorly; pretty-printers are required.

## Verification

Benchmark a tight evaluation loop (10⁸ calls) on line, circle and cubic B-spline. Pass: the erased `Curve` path through `VisitBuiltin` is within 10% of the same loop run through `VisitBuiltin` on a local copy of the geometry (`B0.2/DispatchCopy` beside `B0.2/Dispatch`), and the pure virtual fallback adds at most 5 ns per call over a direct call on the concrete type. The ratio of `VisitBuiltin` to the direct call is reported beside the gate, not gated.

The comparison is with the local copy, revised on 2026-11-01, because that is what isolates the cost of type erasure and dispatch: the same evaluator, called the same way, on the same geometry, held by the wrapper or not. The first form compared with a direct loop on the concrete type, and for the cubic B-spline that ratio read 1.09 on the reference run of 2026-10-28, 1.22 on 2026-10-30 and 1.28 on 2026-11-01. On the run of 2026-11-01 (`B0.2`, [41303397817](https://github.com/onurtuncer/Stereon/actions/runs/41303397817)), `VisitBuiltin` took 16.1 ns, the same loop on a local copy 16.2 ns and a virtual call per evaluation 16.1 ns, against 12.6 ns for the direct loop. On that run, then, the gap did not shrink when the dispatch was made once per loop instead of once per evaluation, nor when the wrapper's storage was replaced by a local copy, so neither is where it comes from as measured. GCC's inlining report on the development build of the twin benchmark (`bench/GeometryBench.cpp`, `-fopt-info-inline-all`, MinGW rather than the runner's Linux build) shows the same decisions on every path: the loop inlined into its caller, the B-spline `Eval` called out of line in all of them, no specialised clone. None of this is a profile of the code the runner executes, so the cause of the gap is not established; code placement such as the hot loop's alignment is a candidate, not a finding (see *Results*).

The bound is absolute because a ratio measures the evaluator rather than the dispatch. A line evaluates in about 2 ns, and in the direct loop the compiler vectorises it, so one indirect call per evaluation costs several times as much whatever the design. For a circle the trigonometry dominates and the ratio is about 1. The benchmark therefore also reports a direct call made one evaluation at a time, for comparison with the fallback.

The same benchmark measures copy cost for buffer sizes 64, 96 and 128 bytes against the actual sizes of the built-in types and fixes the constants.

A tier-1-only test surface wrapped and passed to a boolean must return `NotSupported`; the same surface with `BoundedGeometry` and `ToBSpline` added must take part in the ADR-0001 prototype booleans. The first half is `tests/boolean/TierTests.cpp`: a box whose top face's plane is replaced by the same plane as a tier-1 user type is refused with `NotSupported` by union, cut and common, while the same box with its own plane goes through. The second half waits for the `BoundedGeometry` and `IntersectableGeometry` concepts (ADR-0010).

### Results

`bench/GeometryBench.cpp`, GCC 16.2, -O3 (`docs/geometry.rst`, *Dispatch*), nanoseconds per evaluation:

| Path | Line | Circle |
| --- | --- | --- |
| Direct call on the concrete type | 2.2 | 49 |
| `VisitBuiltin`, loop inside | 2.0 | 46 |
| `Curve`, one virtual call each | 5.6 | 50 |

`VisitBuiltin` is within 10% of direct calls. The fallback adds 3.4 ns for a line and about 1 ns for a circle, both within 5 ns. The same benchmark reports the cubic B-spline row, the scalar direct baseline and the 10⁸-call run.

The reference runs (`B0.2` in `bench/occt/OcctBench.cpp`, the OCCT reference workflow, gcc:16.2 on a hosted ubuntu-24.04 runner; nanoseconds per evaluation, medians of interleaved repetitions):

| Run | Curve | Direct | Local copy | `VisitBuiltin` | Virtual per call | `VisitBuiltin` / copy | `VisitBuiltin` / direct |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2026-10-28 (`8e3aa2e`) | cubic B-spline | 50.9 | — | 55.7 | 55.8 | — | 1.09 |
| 2026-10-30 (`fb7366a`) | cubic B-spline | 22.8 | — | 27.7 | 28.5 | — | 1.22 |
| 2026-11-01 ([41303397817](https://github.com/onurtuncer/Stereon/actions/runs/41303397817)) | line | 0.467 | 0.47 | 0.471 | 0.536 | 1.00 | 1.01 |
| 2026-11-01 | circle | 9.26 | 9.4 | 9.53 | 9.59 | 1.01 | 1.03 |
| 2026-11-01 | cubic B-spline | 12.6 | 16.2 | 16.1 | 16.1 | 0.99 | 1.28 |

Against the revised `VisitBuiltin` clause every curve passes (1.00, 1.01, 0.99). The fallback clause, at most 5 ns per call over the direct call, does not hold for the cubic B-spline on every run: 5.7 ns on 2026-10-30 and 6.4 ns on the run of `main` of 2026-11-01 (41306584204), 4.5 and 3.5 ns on two other runs of that day; it measures the same direct-loop gap as the first form of the `VisitBuiltin` clause did, and is left as written pending the owner's decision. The direct loop's lead on the cubic B-spline, a quarter of an evaluation, is not attributed: the experiments above found no part of it in the dispatch or the wrapper's storage as measured, and the inlining report no difference, but no profile of the runner's code has been taken; that is the next step. It is a cost every loop pays that evaluates a B-spline through `VisitBuiltin` instead of on the concrete type, and it is reported in every reference run.
