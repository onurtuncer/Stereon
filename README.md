# Stereon

[![Linux](https://github.com/onurtuncer/Stereon/actions/workflows/linux.yml/badge.svg?branch=main)](https://github.com/onurtuncer/Stereon/actions/workflows/linux.yml)
[![Windows](https://github.com/onurtuncer/Stereon/actions/workflows/windows.yml/badge.svg?branch=main)](https://github.com/onurtuncer/Stereon/actions/workflows/windows.yml)
[![macOS](https://github.com/onurtuncer/Stereon/actions/workflows/macos.yml/badge.svg?branch=main)](https://github.com/onurtuncer/Stereon/actions/workflows/macos.yml)
[![Lint](https://github.com/onurtuncer/Stereon/actions/workflows/lint.yml/badge.svg?branch=main)](https://github.com/onurtuncer/Stereon/actions/workflows/lint.yml)
[![Python](https://github.com/onurtuncer/Stereon/actions/workflows/python.yml/badge.svg?branch=main)](https://github.com/onurtuncer/Stereon/actions/workflows/python.yml)
[![Benchmarks](https://github.com/onurtuncer/Stereon/actions/workflows/bench.yml/badge.svg?branch=main)](https://github.com/onurtuncer/Stereon/actions/workflows/bench.yml)
[![OCCT reference run](https://github.com/onurtuncer/Stereon/actions/workflows/occt-reference.yml/badge.svg?branch=main)](https://github.com/onurtuncer/Stereon/actions/workflows/occt-reference.yml)
[![Docs](https://github.com/onurtuncer/Stereon/actions/workflows/deploy-docs.yml/badge.svg?branch=main)](https://onurtuncer.github.io/Stereon/)
[![codecov](https://codecov.io/gh/onurtuncer/Stereon/branch/main/graph/badge.svg)](https://codecov.io/gh/onurtuncer/Stereon)
[![License: MPL 2.0](https://img.shields.io/badge/License-MPL_2.0-brightgreen.svg)](LICENSE)
[![C++26](https://img.shields.io/badge/C%2B%2B-26-blue.svg)](docs/adr/0016-cpp26-adoption-and-contracts.md)

**A modern C++26 boundary-representation (B-rep) geometric modeling kernel.**

Stereon is an open-source solid modeling kernel: NURBS and analytic geometry, B-rep topology, robust Boolean operations, STEP data exchange, watertight meshing for analysis and certified queries for machining. It is built from scratch on C++26 with value semantics, type-erased geometry and thread safety by construction.

> **Status: pre-alpha, version 0.0.4.** Stereon is finishing Phase 0 (foundations), and the v0.1 geometry gates are implemented. APIs change without notice, and nothing here is ready for production use. See the [roadmap](#roadmap) for what exists today and what comes next.

---

## The name

*Stereon* (στερεόν) is the ancient Greek word for "a solid": the neuter of *stereós*, "firm, solid". It is the word Euclid uses in the *Elements*, where Books XI–XIII treat solid geometry. Book XI opens with two definitions:

> Στερεόν ἐστι τὸ μῆκος καὶ πλάτος καὶ βάθος ἔχον. Στερεοῦ δὲ πέρας ἐπιφάνεια.
>
> *A solid is that which has length, breadth and depth. And the boundary of a solid is a surface.*
>
> — Euclid, *Elements* XI, Defs. 1–2 (ed. J. L. Heiberg, *Euclidis Elementa*, vol. IV, Leipzig: Teubner, 1885)

The second definition is the idea behind boundary representation: a solid is fully described by the surfaces that bound it. That is how Stereon models shapes. The same root survives in *stereometry*, the old name for the measurement of solids.

## Why Stereon

Mature kernels such as OpenCASCADE carry decades of design history: reference-counted handle graphs, global state, and tolerance behaviour that is hard to predict. Stereon aims to be narrower and sharper:

- **Robust core operations.** Booleans on NURBS and analytic B-rep with a documented tolerance model, exact predicates for every topological decision, and a public validity checker. An operation either returns a valid solid or a structured error — never a silently broken shape.
- **A modern, embeddable API.** Value types, no global state, `std::expected` error handling, C++20 modules, and Python bindings.
- **Parallel and deterministic.** Operations run on a `std::execution` task graph and give bit-identical results regardless of thread count.
- **Built for engineering analysis.** Watertight surface meshes for CFD/FEA, and (planned) shape sensitivities for design optimisation.
- **Built for machining (planned).** B-rep to tri-dexel conversion for in-process stock, ray, distance, contact and clearance queries, and tool offsets, each with certified error bounds and explicit unresolved cases ([ADR-0044](docs/adr/0044-brep-to-tri-dexel-conversion.md)). Toolpath strategy and postprocessing stay in applications.

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

Third-party dependencies are fetched by CMake at pinned versions (see [`cmake/StereonDependencies.cmake`](cmake/StereonDependencies.cmake)): [stdexec](https://github.com/NVIDIA/stdexec) (until the standard library ships `std::execution`) and [xsimd](https://github.com/xtensor-stack/xsimd) (until it ships `std::simd`) always; Catch2 and RapidCheck with the tests; Google Benchmark with `STEREON_BUILD_BENCH`; nanobind with `STEREON_BUILD_PYTHON`; [CppAD](https://github.com/coin-or/CppAD), the AD backend, with `STEREON_WITH_CPPAD` for its integration tests only (EPL-2.0, never linked into a library; [ADR-0045](docs/adr/0045-cppad-ad-backend.md)).

### Build and test

```bash
git clone https://github.com/onurtuncer/Stereon.git
cd Stereon
cmake --preset gcc16-release
cmake --build --preset gcc16-release
ctest --preset gcc16-release
```

Other presets: `gcc16-debug` (contracts enforced, checked handles), `gcc16-observe` (contracts logged, for corpus runs), `gcc16-asan` and `gcc16-tsan` (Linux and macOS), `gcc16-windows` (static GCC runtime) and `bench`. Release presets keep API-boundary contracts and compile out hot-path ones (ADR-0016).

Add `-DSTEREON_BUILD_PYTHON=ON` to any preset to build the Python module (`import stereon`).

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
├── tools/           # source checks run by CI (layering, epsilon ban)
└── docs/            # Sphinx manual (placeholders), built with -DSTEREON_BUILD_DOCS=ON
    └── adr/         # architecture decision records
```

## Roadmap

Each phase ends with an exit gate checked in CI. The [detailed roadmap and versioning plan](docs/planning/roadmap-and-versioning.md) embeds modelling subgates, acceptance evidence and compatibility rules into this release sequence.

| Release | Scope and milestone gates | Status |
| --- | --- | --- |
| Phase 0 / 0.0.x | Foundations (M1): build/CI, math, predicates, type erasure, controlled-edit prototype, debug viewer, diagnostics and minimal persistence prototype | 🚧 Exit gate nearly met |
| v0.1 "Geometry" | Analytic/NURBS evaluation, fitting and projection; geometry persistence and numerical-domain checks (M1/M4) | 🚧 Implemented, awaiting the Phase 0 gate |
| v0.2 "Intersections" | Curve/curve, curve/surface and surface/surface services; explicit tangency, overlap, coincidence and unresolved outcomes (M2/M4) | Planned |
| v0.3 "Solids" | Planar construction, sectioning and imprinting (M2); curved/trimmed construction, primitives, extrude/revolve/sweep/loft, sewing, history and native persistence (M4/M6); certified discretisation, ray–face queries and tri-dexel conversion of analytic solids (M7) | Planned |
| v0.4 "Planar Booleans" | Validated union/cut/common on the declared planar/polyhedral subset, with lineage and tags (M3) | Planned |
| v0.5 "Curved Booleans" | Bounded analytic/trimmed-spline Boolean integration; seams, poles, tangencies and certificates (M4); tri-dexel conversion and closest-point queries on trimmed B-spline faces (M7) | Planned |
| v0.6 "Booleans" | Broad declared real-part corpus, robustness, comparison, tolerance-growth and determinism gates (M3/M4) | Planned |
| v0.7 "Meshing" | Derived display/analysis meshes, tagged engineering workflows, native migrations and usable Python operation calls (M4/M6); tri-dexel fields with interval operations, tool-contact and clearance queries, B5.7 gate of ≥ 3× OCCT (M7) | Planned |
| v0.8 "Exchange" | STEP AP242 with healing, derived watertight meshes, native migrations and usable Python/tagged engineering workflows (M4/M6) | Planned |
| v0.9 "Finishing" | Bounded fillets, chamfers, shelling/thickening and shape sensitivities; history/tag integration (M5/M6); offset surfaces for tool-centre geometry (M7) | Planned |
| v1.0 "Stable kernel" | Stable public API, compatibility policy and supported-domain release qualification across M1–M7 | Planned |

**Where things stand (0.0.4).** Phase 0's exit gate has two lines left: accepting [ADR-0007](docs/adr/0007-threading-and-determinism.md), whose verification landed with the Boolean prototype's M10, and the OCCT differential tests of the prototype's Booleans (M11). Every blocking CI job is green on Linux, Windows and macOS. The Boolean prototype (union, cut and common of boxes and cylinders, analytic intersections with interval certificates, classification by rays and a bounded winding number, operation history) runs through M10.

For v0.1, every gate in the [acceptance ledger](docs/planning/v0.1-acceptance-ledger.md) has its implementation: analytic and NURBS evaluation with derivatives, global projection, exact NURBS conversion, curve and surface fitting, trimmed curves and surfaces, extrusion and revolution surfaces, and `.stn` persistence of the geometry with its numerical policy. The B-spline operations OCCT users expect are also in: knot insertion and removal, degree elevation, clamping, periodic form, concatenation, continuity queries, and SIMD batch and grid evaluation that return the point evaluator's results bit for bit. The first derivative products are in too: forward and reverse products (JVP, VJP) of B-spline curve points with respect to poles and weights (ADR-0021, stage 2), and every public operation of the curve, surface, intersection, topology, construction and Boolean libraries declares its derivative capability. Every benchmarked row is at or above parity with OCCT on the reference run; the newest rows (batch, grid and degree elevation) wait for their first reference run. The v0.1 tag follows the Phase 0 gate, together with the release report the roadmap requires.

v0.2 carries the Phase 2 intersection gate. Phases 4–6 span several minor releases; phase numbers are not product version numbers. Milestones M1–M7 are dependency checkpoints (M7, machining services, is Stereon's own; M1–M6 follow Stroud), while the [capability matrix's R0–R5](docs/planning/stereon-capability-matrix.md) track broader long-term coverage. Drafting, variable-radius blends and full parametric history remain post-v1.0. Product, native-format and C ABI versions are separate; a planning change does not bump the current kernel version.

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


---

## Community

- **Bug Tracker**: [GitHub Issues](https://github.com/onurtuncer/Stereon/issues)
- **Discussions**: [GitHub Discussions](https://github.com/onurtuncer/Stereon/discussions/)

---

## Author

**Prof. Dr. Onur Tuncer**  
Aerospace Engineer, Researcher & C++ Systems Developer  
Istanbul Technical University · [onur.tuncer@itu.edu.tr](mailto:onur.tuncer@itu.edu.tr)

<p align="left">
  <img src="assets/itu_logo.png" width="180" alt="Istanbul Technical University"/>
</p>
