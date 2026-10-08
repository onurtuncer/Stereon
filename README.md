# Stereon

**A modern C++26 boundary-representation (B-rep) geometric modeling kernel.**

Stereon is an open-source solid modeling kernel: NURBS and analytic geometry, B-rep topology, robust Boolean operations, STEP data exchange and watertight meshing for analysis. It is built from scratch on C++26 with value semantics, type-erased geometry and thread safety by construction.

> **Status: pre-alpha.** Stereon is in Phase 0 (foundations). APIs change without notice, and nothing here is ready for production use. See the [roadmap](#roadmap) for what exists today and what comes next.

---

## Why Stereon

Mature kernels such as OpenCASCADE carry decades of design history: reference-counted handle graphs, global state, and tolerance behaviour that is hard to predict. Stereon aims to be narrower and sharper:

- **Robust core operations.** Booleans on NURBS and analytic B-rep with a documented tolerance model, exact predicates for every topological decision, and a public validity checker. An operation either returns a valid solid or a structured error — never a silently broken shape.
- **A modern, embeddable API.** Value types, no global state, `std::expected` error handling, C++20 modules, and Python bindings.
- **Parallel and deterministic.** Operations run on a `std::execution` task graph and give bit-identical results regardless of thread count.
- **Built for engineering analysis.** Watertight surface meshes for CFD/FEA, and (planned) shape sensitivities for design optimisation.

Out of scope for v1.0: parametric feature history, sheet metal, drafting, variable-radius blends.

## Design at a glance

| Principle | How |
| --- | --- |
| External polymorphism | Geometry types (`Line`, `BSplineCurve`, `Cylinder`, `BSplineSurface`, …) are plain structs with no base class. `Curve` and `Surface` wrappers erase them; any type that satisfies the concept plugs in. A closed `std::variant` fast path keeps hot loops free of virtual dispatch. |
| Data-oriented topology | Vertices, edges, coedges, loops, faces, shells and solids live in structure-of-arrays stores addressed by typed 32-bit handles. Radial coedge rings support non-manifold intermediate results. |
| Immutable shapes | A mutable `ShapeBuilder` freezes into an immutable `Shape`. Copy-on-write pages make snapshots and undo cheap and make shapes safe to share across threads. |
| Explicit tolerances | A global model resolution plus per-edge and per-vertex tolerances that only grow through explicit, logged operations. |
| History maps | Every operation reports which input entities generated, modified or deleted which output entities. |
| C++26 throughout | Contracts on every geometric operation, static reflection for serialisation and bindings, `std::execution` for parallelism, `std::simd` in evaluators. |

Design decisions are recorded as ADRs in [`docs/adr/`](docs/adr/).

## A taste of the API

*Target API for v0.6; not yet implemented.*

```cpp
import Stereon;
using namespace Stereon;

int main() {
    Context ctx;                                        // no global state

    auto box  = MakeBox(ctx, {0, 0, 0}, {40, 30, 20});
    auto hole = MakeCylinder(ctx, Axis{{20, 15, -1}, {0, 0, 1}}, 6.0, 22.0);

    std::expected<Shape, KernelError> part = BooleanCut(ctx, *box, *hole);
    if (!part) {
        std::println(stderr, "cut failed: {}", part.error().GetMessage());
        return 1;
    }

    std::println("volume = {:.3f} mm^3", ComputeMassProperties(*part).Volume);
    WriteStep(*part, "bracket.step");
}
```

## Building

### Requirements

- **GCC 16 or newer** (reference compiler). As of 2026 it is the only compiler with the C++26 reflection and contracts support Stereon relies on. Clang and MSVC run in CI as non-blocking until they catch up.
- CMake 3.30+ and Ninja
- Python 3.11+ (optional, for bindings and test tooling)

Third-party dependencies are fetched by CMake: [stdexec](https://github.com/NVIDIA/stdexec) (until the standard library ships `std::execution`), [xsimd](https://github.com/xtensor-stack/xsimd) (fallback for `std::simd`), Catch2, RapidCheck, Google Benchmark and nanobind.

### Build and test

```bash
git clone https://github.com/onurtuncer/Stereon.git
cd Stereon
cmake --preset gcc16-release
cmake --build --preset gcc16-release
ctest --preset gcc16-release
```

Other presets: `gcc16-debug` (contracts enforced, checked handles), `gcc16-observe` (contracts logged, for corpus runs), `gcc16-asan` and `gcc16-tsan` (Linux and macOS), `gcc16-windows` (static GCC runtime) and `bench`.

On Windows, build from an MSYS2 UCRT64 shell with `mingw-w64-ucrt-x86_64-gcc`, `-cmake` and `-ninja` installed.

## Repository layout

```
stereon/
├── src/
│   ├── core/        # math types, intervals, boxes, handles, arenas, errors
│   ├── robust/      # exact predicates, interval and expansion arithmetic
│   ├── curve/       # 2D/3D analytic and NURBS curves
│   ├── surface/     # analytic, NURBS and procedural surfaces
│   ├── intersect/   # curve-curve, curve-surface, surface-surface
│   ├── topo/        # B-rep store, Euler operators, validity checker
│   ├── build/       # primitives, extrude, revolve, sweep, loft, sewing
│   ├── boolean/     # union, cut, common
│   ├── blend/       # fillets, chamfers, shelling
│   ├── mesh/        # display and analysis-grade tessellation
│   └── io/          # STEP, STL, OBJ, glTF, native format
├── python/          # nanobind bindings
├── viewer/          # debug viewer for any intermediate object
├── tests/           # unit, property-based, fuzz and regression tests
├── corpus/          # geometry test corpus (regression cases never deleted)
├── bench/           # benchmark suite
└── docs/
    └── adr/         # architecture decision records
```

## Roadmap

Each phase ends with an exit gate that is checked in CI.

| Release | Scope | Status |
| --- | --- | --- |
| Phase 0 | Build system, CI, core math, robust predicates, type-erasure scaffolding, debug viewer | 🚧 In progress |
| v0.1 "Geometry" | NURBS and analytic curves/surfaces: evaluation, fitting, projection | Planned |
| v0.3 "Solids" | B-rep topology, validity checker, primitives, extrude, revolve, sweep, loft | Planned |
| v0.6 "Booleans" | Robust union, cut and common on real parts | Planned |
| v0.8 "Exchange" | STEP AP242 read/write with healing, watertight meshing, Python bindings | Planned |
| v1.0 | Fillets, chamfers, shelling, shape sensitivities, stable API | Planned |

## Testing philosophy

A kernel is only as good as the ugly cases it has survived. Stereon relies on:

- **Differential testing** against OpenCASCADE and Manifold, plus analytic volumes and areas.
- **Property-based tests**, e.g. `vol(A ∪ B) + vol(A ∩ B) = vol(A) + vol(B)`.
- **Fuzzing** of the STEP parser and of random primitive placements.
- **A growing failure corpus.** Every bug becomes a regression case that is never deleted.
- **Determinism checks.** The same input must give bit-identical output on any number of threads.

## Contributing

Contributions are welcome once the Phase 0 foundations settle. Before opening a pull request, please read:

- [`CONTRIBUTING.md`](CONTRIBUTING.md): coding standards, commit conventions, and the contributor agreement
- [`docs/adr/`](docs/adr/): accepted design decisions. A PR that contradicts an accepted ADR needs a new ADR.

Every change must include a test that fails before it and passes after it. AI-assisted contributions are welcome under the same rule; see [ADR-0015](docs/adr/0015-ai-contribution-policy.md) for the policy.

## License

Stereon is licensed under the [Mozilla Public License 2.0](LICENSE).

## Citing

If you use Stereon in academic work, please cite the repository until a paper is available:

```bibtex
@software{stereon,
  title  = {Stereon: a modern C++26 B-rep geometric modeling kernel},
  author = {Tuncer, Onur and contributors},
  year   = {2026},
  url    = {https://github.com/onurtuncer/Stereon}
}
```

## Acknowledgements

Stereon builds on ideas from decades of solid modeling research, in particular *The NURBS Book* (Piegl & Tiller), *Shape Interrogation for Computer Aided Design and Manufacturing* (Patrikalakis & Maekawa), *Boundary Representation Modelling Techniques* (Stroud), and Shewchuk's robust geometric predicates.
