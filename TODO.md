# Phase 0 — Foundations: TODO

See the [roadmap and versioning plan](docs/planning/roadmap-and-versioning.md) for the M1 foundation milestone, early persistence staging and downstream release gates. This checklist remains the Phase 0 implementation gate.

**Goal:** a building, tested skeleton of Stereon with robust predicates, type-erased geometry scaffolding and the Phase 0 ADRs settled.
**Duration:** ~12 weeks (8–10 with AI-assisted development).
**Exit gate:** see the [last section](#exit-gate).

---

## 1. Repository and governance (week 1)

- [x] Create the `onurtuncer/Stereon` repository (public, kept under the personal account)
- [x] Add `LICENSE` (MPL-2.0, per ADR-0008) and the MPL header template
- [x] Add `README.md`, `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, `SECURITY.md`
- [x] Add `docs/adr/` with ADR-0000 template and ADRs 0001–0016
- [x] PR template with test checklist and AI-disclosure checkbox (ADR-0015)
- [x] Issue templates: bug (with native-file or journal attachment), feature, corpus case
- [x] `CODEOWNERS` with module owners for Tier B/C libraries (ADR-0015)
- [x] Branch protection on `main`: required CI checks, one approving review, linear history (admins may bypass while there is a single maintainer)
- [x] Set the GitHub About description and topics
- [ ] Choose a logo direction from the four concepts; add the mark and wordmark as SVG under `docs/brand/` and use it in the README header

## 2. Build system and toolchain (weeks 1–2)

- [x] Top-level `CMakeLists.txt` (CMake ≥ 3.30, Ninja), C++26, modules enabled
- [x] One CMake target per library (`stereon::core`, `stereon::robust`, …) with private include visibility. All eleven libraries exist; only `src/<library>/include/` is public
- [x] `CMakePresets.json`: `gcc16-debug`, `gcc16-release`, `gcc16-observe`, `gcc16-asan`, `gcc16-tsan`, `gcc16-windows`, `bench`
- [x] Contract evaluation semantics per preset: enforce (debug), observe (CI corpus), ignore/enforce-boundary (release) (ADR-0016). Release enforces `STEREON_PRE`/`POST`/`ASSERT` and compiles out the `STEREON_HOT_*` contracts
- [x] Dependency fetching (CPM or FetchContent): stdexec, xsimd, Catch2, RapidCheck, Google Benchmark, nanobind
- [x] Feature-test-macro wrappers: `Stereon::Exec` (std::execution ↔ stdexec), `Stereon::Simd` (std::simd ↔ xsimd). `Exec` is a header, not a module, because GCC 16.2 and MSVC crash on stdexec in modules
- [x] Include-graph check enforcing library layering (fails CI on an upward dependency): `tools/check_layering.py`
- [x] `.clang-format` and `.clang-tidy`, including a rule banning literal epsilon comparisons in `topo/`, `boolean/`, `mesh/` (ADR-0001, ADR-0002). The ban is `tools/check_epsilon.py`, since clang-tidy cannot parse the GCC build; CI checks formatting with clang-format 22.1.8
- [x] Licence-header check (REUSE or equivalent) (ADR-0008)
- [x] Function complexity and length limits per tier (ADR-0022): `tools/check_complexity.py` in the lint job and CTest, with the self-test the ADR asks for; `tools/requirements.txt` pins its lizard and CI's clang-format
- [ ] Remove `stereon_anchor_contracts()` and `ContractAnchor()`, restore `RelWithDebInfo` for `gcc16-observe` and hot-path contracts for `gcc16-tsan` once the known-failing `stereon_contract_probe` test turns green, i.e. GCC checks contracts imported from modules (docs/toolchain/gcc16-feature-probe.md)

## 3. Continuous integration (week 2)

- [x] Linux GCC 16: debug + release, blocking
- [x] Linux GCC 16: ASan/UBSan and TSan jobs, blocking
- [x] Windows GCC 16 via MSYS2 UCRT64 (`msys2/setup-msys2`), blocking
- [x] macOS GCC 16 (Homebrew), blocking
- [x] Clang (latest) and MSVC (latest): non-blocking, reported
- [x] Coverage report (gcov/llvm-cov) published per PR
- [x] Benchmark job on `main` with results stored for trend tracking. The `trend` job of `.github/workflows/bench.yml`: `stereon_bench` with five repetitions on every push to main, one run file per commit on the `bench-data` branch (`tools/bench_trend.py`), regressions over 10 % against the run before as warnings and in the job summary, and the trend page published with the manual under `benchmarks/trend/`
- [x] Nightly job: full test suite + determinism check (1, 2, 8, 32 threads) (ADR-0007). `.github/workflows/nightly.yml`: the suite on gcc16-debug and gcc16-release, and `stereon_determinism` per configuration in its own process, twice, compared byte for byte (`docs/testing.rst`, *Determinism*)
- [x] Dependabot keeps the GitHub Actions and the pip pins (`docs/requirements.txt`, `tools/requirements.txt`) current: monthly, grouped
- [ ] Check the CMake FetchContent pins (stdexec, xsimd, Catch2, RapidCheck, Google Benchmark, nanobind) and the `gcc:16.2` image against upstream at each phase gate; Dependabot cannot read them

## 4. `stereon::core` (weeks 2–5)

- [x] Fixed-size `Vec2/3`, `Point2/3`, `Mat3/4`, `Transform` with constexpr operations. They are the `double` instances of scalar-generic `Basic*` templates over the `Scalar` concept (ADR-0021), checked with a forward-mode dual number in tests
- [x] `Interval` with outward rounding. Bounds step one representable value outward on the bit pattern, so the arithmetic is constexpr and needs no rounding-mode state
- [x] `Box2`, `Box3` (bounding boxes) with union, intersection and interval-based tests. A box is one `Interval` per axis, so every test is the per-axis interval test
- [x] Typed handles `Id<Tag>` (32-bit), one distinct type per entity kind, no implicit integer conversions (ADR-0004). `Id<EntityKind>` with the kind as the tag, so diagnostics and the C ABI can name it
- [x] Generational builder handles (index + generation) with free lists (ADR-0004). `BuilderId` and `GenerationTable`; slots released in a session become reusable at `BeginSession()`, which keeps indices stable across one edit
- [x] `STEREON_CHECKED_HANDLES` mode (on in debug, asan, fuzzing presets): handles carry store ID + generation; every access checks store, bounds and generation (ADR-0004 safeguards). `GenerationTable::Check()`; the report and contract violation come from `Handle.cpp`, a module implementation unit, so they fire in every importer
- [x] `static_assert` + CI check that release-build handles are exactly 4 bytes (`Handle.cppm`; the release and Windows CI presets compile it)
- [x] Named accessors for packed orientation bits (`IsReversed(h)`, `BaseOf(h)`); no manual masking. `Oriented<Kind>` keeps the sense in the top bit and `Reversed(h)` flips it
- [x] Paged arena / copy-on-write page container (`shared_ptr<const Page>`, 256 entries) (ADR-0005). `PagedArena<T, N>`: a copy shares every page, the first write to a page copies it unless the arena holds it alone, and hot loops read one page span at a time; `bench/ArenaBench.cpp` runs the 128/256/1 024 sweep on a 10⁶-record table
- [x] `KernelError`, `ErrorCode`, `std::expected` aliases, `STEREON_TRY` macro (ADR-0006). `KernelError` carries code, message and tolerances; entities, `ValidityReport` and journal ID follow with handles, the checker and the journal
- [x] `Context`: resolution, tolerance ceiling, scheduler, logger, journal sink (ADR-0001, ADR-0007). A value with `With...` modifiers; also the size box, the tolerance-growth policy (ADR-0006), the validity level (ADR-0018) and the compaction threshold (ADR-0004). The scheduler is a kind and a thread count, since stdexec types cannot appear in module code (docs/toolchain/gcc16-feature-probe.md)
- [x] Logging and operation journal (inputs + parameters, replayable). Kernel-owned `LogSink` and `JournalSink` (ADR-0023), `MemoryJournal`, `JournalId` on `KernelError`; replay of a recorded entry lands with the first operation (the ADR-0001 prototype in section 10)
- [x] Unit and property tests for all of the above. Every core item has Catch2 cases and RapidCheck properties in `tests/core/`
- [x] Sensitivity foundation types (ADR-0021 stage 1). Names stay provisional until ADR-0017 is applied at the first slice; no AD backend in any of them:
  - [x] Parameter identities: a typed `ParameterId` with its kind and unit (lengths in mm, angles in radians, per ADR-0001) and scaling convention, plus the ordering convention for parameter vectors and seeds (§3, §7). `ParameterId`, `ParameterSet` with per-parameter reference scales in `Stereon.Core:Sensitivity`
  - [x] Derivative statuses `Valid`, `NearEvent`, `Nondifferentiable`, `NotSupported`, `Indeterminate`, and a `SensitivityResult` carrying the status, the branch/topology signature, affected entities, tolerance policy, residuals and conditioning diagnostics; returned as `std::expected<SensitivityResult, KernelError>`, with no usable payload for `NotSupported`, `Indeterminate` and `Nondifferentiable` (§5). `GetValues()` is empty, never zeros, for the three
  - [x] Capability declarations: a per-operation descriptor of parameter dependencies, correspondence rules, differentiability assumptions and supported products (JVP, VJP, dense Jacobian), queryable before evaluation; `NotSupported` by default, never a silent zero (§5, Decision). `OperationCapability` and one catalogue per library (`src/<library>/Capabilities.cppm`); every product NotSupported today
  - [x] Correspondence types: a tracked sample naming its source entity, orientation, local coordinates and tracking convention (fixed surface coordinates, fixed edge parameter, moving trim coordinates, unique projection), independent of entity indices so `Compact()` cannot redefine it (§4). Entity references use the §7 handles once they exist. `TrackedSample` names the snapshot's store and the slot together, so a renumbered snapshot never matches it
  - [x] Sensitivity context: caller-owned, holds seeds, adjoints and solver workspaces for one evaluation, separate from `Context`, never stored in a `Shape`; two contexts share no mutable state (§3)
  - [x] Registration check in tests: adding a geometric operation without a capability declaration fails the build or the test suite (ADR-0021 verification). `tools/check_capabilities.py`, in lint and CTest

## 5. `stereon::robust` (weeks 3–6)

- [x] Adaptive-precision `Orient2D`, `Orient3D`, `InCircle`, `InSphere` (Shewchuk-style). `Stereon.Robust`: stage A with Shewchuk's bounds, stage B on the rounded differences (exact when their tails vanish), stage C for the orientations, then the exact cofactor form on the raw coordinates; results are `Sign`, never doubles
- [x] Expansion arithmetic primitives for custom predicates. `TwoSum`, `TwoDiff`, `TwoProduct`, `Split`; `Expansion<N>` with capacity-typed `Sum`, `Scaled`, `Product` and operators, all constexpr; the library is built with `-ffp-contract=off` (PUBLIC)
- [x] Interval-filtered predicate helpers with exact fallback. `SignOf(Interval)` and `FilteredSign(estimate, exact)`
- [x] Adversarial test suite (Shewchuk cases + 10⁶ random near-degenerate cases). The Kettner grid, an `__int128` oracle on integer and dyadic inputs, exactly degenerate integer constructions (collinear, coplanar, cocircular, cospherical) with one-unit perturbations; 20 000 cases per predicate in CI, the 10⁶ run is the hidden test `stereon_robust_tests "[exhaustive]"`
- [x] Benchmarks: filtered vs naive floating point (`bench/RobustBench.cpp`, B0.1): random, nearly degenerate and exact-stage costs per predicate

## 6. Geometry scaffolding (weeks 4–8)

- [x] Concepts `CurveGeometry`, `Curve2dGeometry`, `SurfaceGeometry` (ADR-0003). Tier 1; tiers 2 and 3 arrive with the intersection work
- [x] Type-erased `Curve`, `Curve2d`, `Surface` with 96-byte small-buffer storage and shared payloads. `Stereon.Core:SmallBuffer`; types too large for the buffer are shared immutable
- [x] Two toy types per concept: `Line` + `Circle`, `Plane` + `Cylinder`. Also `Line2d` + `Circle2d` for `Curve2dGeometry`; all scalar-generic (ADR-0021)
- [x] `VisitBuiltin` fast path generated with expansion statements and pack indexing (ADR-0016). Within 10% of direct calls in `bench/GeometryBench.cpp`; the virtual fallback is 2.6x direct for a line, 1.03x for a circle
- [x] Type-pair dispatch table for binary operations (stub intersection handlers). `Stereon.Intersect`: compile-time tables over `BuiltinCurveTypes` and `BuiltinSurfaceTypes`, `NotImplemented` stubs, `NotSupported` for user types
- [x] Contracts on all evaluator entry points (domain, finiteness). Plain branches reporting from `EvaluatorChecks.cpp`, because GCC 16.2 drops contracts in imported templates; domain checks honour `Period`/`Periods`; hot-path, so debug and sanitizer presets only
- [x] Reflection-based serialiser prototype for geometry structs (ADR-0012 groundwork). `Serialized`, `Deserialized<T>` (ParseError on malformed input) and `ToJson` in `Stereon.Core:Serialize`
- [x] Property tests comparing fast path and erased path on every type. Bitwise agreement for every operation of every built-in curve and surface

## 7. Topology prototype (weeks 6–10)

- [x] SoA stores for vertex, edge, coedge (with radial ring), loop, face, shell, solid, following Golovanov's logical entity model (ADR-0004). `Stereon.Topo:Records`: one `PagedArena` record table per kind with the slot state in the record; splitting records into per-field arrays is left to the ADR-0004 benchmark
- [x] Edge record per ADR-0009: curve kind (`Intersection`, `Exact`, `Trim`, `Boundary`, `Spatial`, `Degenerate`), tolerance, stored primary coedge, own or cached 3D curve slot with certificate; pcurves on coedges on a shared parameter. `ShapeBuilder::CacheEdgeCurve` builds the cache from the primary side; its certificate is `Uncertified` until interval evaluation of S(c(t)) exists (ADR-0009, section 4), and `Certified` is reserved for it; rule S17 checks the cache fields
- [x] Cache invalidation of the 3D edge curve when a pcurve changes (stub is enough for Phase 0). `AppendCoedge`, the removals, `SetPrimaryCoedge` and `SetSurface` drop the caches of the edges they touch; `Commit` warns with `CacheDropped`
- [x] `ShapeBuilder` → `Freeze()` → immutable `Shape` with stable indices and alive bitsets; separate `Compact()` returning the old→new remap table (ADR-0004). `Freeze()` certifies a `TopologySnapshot` (ADR-0038) and keeps every index; the `Remove` operations leave tombstones and keep every link consistent; `ShapeBuilder::Open(shape)` edits an existing `Shape` sharing its pages; `Compact()` renumbers and returns the `RemapTable`
- [x] Adjacency views: `FacesOf(edge)`, `EdgesOf(vertex)`, loop walks. `AllEntities`, `LoopsOf`, `CoedgesOf`, `RadialCoedgesOf`, `FacesOf(shell)` and `ShellsOf` ranges (ADR-0036); `EdgesOf(vertex)` from the maintained vertex-incidence index (rule S16 checks it), `FacesOf(edge)` as the eager distinct result
- [x] `Structural`-level validity checker with rule IDs and a `ValidityReport`, naming offending entities with `Describe()` (ADR-0018); `Geometric` level once edge certificates exist. Rules S01–S17, bounded walks over the raw tables; the `Geometric` level is still to come
- [x] `OperationResult` (shape, history, tolerance changes, warnings) as the return type of shape-producing operations (ADR-0006). `Commit(builder, ctx)` and `CompactShape(shape, ctx)`: a History of Unchanged, Modified, Generated and Deleted relations (split and merge events join with the Euler operators), `SetTolerance` changes coalesced per entity under the Context's ceiling and growth policy, `ToleranceGrown` and `CacheDropped` warnings, and the level checked
- [x] `Describe(shape, handle)` for every entity kind, used in logs, contract messages and checker reports
- [x] Debug builds run the validity checker after every public operation (ADR-0004 safeguards). At publication (`Freeze`, `Commit`, `Compact`, `CompactShape`) always; checked builds also validate the input of `Open`, `Commit` and `CompactShape` again, and every operation checks its History before publishing. A builder's single steps are not checked, since a valid edit passes through states that are not (`docs/topology.rst`, *Operation results*). `Geometric` follows with the ADR-0018 level
- [x] Build a box (`Exact` edges) and a cylinder (`Exact` + `Boundary`/seam edges) by hand as test bodies. `tests/topo/Bodies.hpp`, with prisms, sheets, a holed face, a non-manifold fan, a wire and a degenerate edge

## 8. Debug viewer v0 (weeks 5–9)

Per [ADR-0035](docs/adr/0035-debug-viewer.md): objects reach the viewer as `.stn` dump files; the viewer is GLFW + OpenGL 3.3 core + Dear ImGui.

- [x] `Stereon.Mesh:Sampling`: curves to polylines and surfaces to iso-lines within a display tolerance (chord and angle) kept apart from modelling tolerances; NURBS sampled through every knot; unbounded geometry over a given range. `SampleCurve`, `SampleSurface` and `Breakpoints` in `src/mesh/Sampling.cppm`; the chord is checked at seven interior points per segment, so the bound is for display, not certified
- [x] `Stereon.Io:Dump`: one `.stn` file per call in `$STEREON_DUMP_DIR` (nothing written when unset), the geometry chunks plus an optional `DBUG` chunk (name, source location, points, polylines, boxes, triangle soups, labels); written to a temporary file and renamed. `Stereon.Io:Debug` holds the chunk's layout; a `Shape` is dumped through `ShapeArchive` (geometry plus overlay primitives) until the format has topology chunks
- [x] `STEREON_DUMP(name, obj)` in `<Stereon/Io/Dump.hpp>`, and exported non-template `Dump(name, obj)` overloads callable from a GDB or LLDB breakpoint. `Dump(const char*, ...)` for `Curve`, `Curve2d`, `Surface` and `Shape`
- [x] `viewer/scene`: headless library turning dump files into display primitives, with Catch2 tests that open no window. Module `Stereon.Viewer.Scene`; `tests/viewer/SceneTests.cpp` runs in every CI job
- [x] `viewer/app` (`stn-view`) behind `STEREON_BUILD_VIEWER` (off by default): orbit/pan/zoom, fit all, object list with visibility and source location, live reload of the dump directory. Plus `--screenshot` for the manual's figures; OpenGL through a hand-written loader (`viewer/app/Gl.hpp`), so no generated loader is a dependency
- [x] Geometry overlays: parameter ticks, direction arrows, control polygons and nets, knots; a 2D panel for `Curve2d` and pcurves in parameter space. The 2D panel draws the plane curves of the scene; a pcurve's surface domain needs the topology chunks
- [ ] Topology overlays once §7 lands: edges coloured by curve kind, coedge orientation, vertex and edge tolerances as spheres and tubes, `Describe()` labels, `ValidityReport` highlighting by rule ID. Partly: edges coloured by kind, vertex tolerances as boxes and `Describe()` labels come with a `Shape` dump; coedge orientation and `ValidityReport` highlighting are open
- [x] GLFW, Dear ImGui and the GL loader declared with licences and admitted by `tools/check_dependencies.py` only under `viewer/`; a check that nothing under `src/` references them (ADR-0034); CI builds the viewer and runs the `viewer/scene` tests. The Linux viewer job also renders the manual's screenshots under Xvfb
- [x] Debugger pretty-printers for GDB (Python) and LLDB: handles resolved to the entity (e.g. `Edge#4711 [Line, tol 1e-5, faces 12|37]`), `Interval`, `Box3`, `Curve`, `Surface`. Done except resolving a handle to its entity, which needs the topology stores: `tools/debugger/stereon_gdb.py` and `stereon_lldb.py` print vectors, intervals, boxes, handles, the erased wrappers (built-in and user types), arenas, journal ids, signs and expansions as `std::format` does. Handles of a live shape now show their entity in checked builds: snapshots register their store, and the printers call `stereon_describe_handle()` (`tests/debugger/check_printers.py --checked-handles`)
- [x] CI test for the pretty-printers using scripted GDB and LLDB sessions. `tests/debugger/check_printers.py` against one expected table for both debuggers; the Linux CI jobs install gdb and lldb
- [x] `stn-dump` prototype: JSON view of any `.stn` file, dumps included (geometry now, `Shape` once it exists), for diffing two versions (ADR-0012 groundwork). `tools/stn-dump`, with `--pretty` for diffs; `Shape` dumps show their overlays, not yet their topology records

## 9. Testing infrastructure (weeks 2–10)

- [x] Catch2 + RapidCheck wiring; test naming and tagging conventions. One executable per library, cases tagged `[<library>]` plus topic tags and `[property]` for `rc::prop` tests ([`tests/CMakeLists.txt`](tests/CMakeLists.txt))
- [x] libFuzzer harness skeleton (for later STEP and primitive-placement fuzzing). `tests/fuzz/`: the `.stn` reader and the checked B-spline factory, behind libFuzzer's entry point; a standalone driver for GCC (seed corpus and seeded mutations, contract violations counted), a CTest smoke pass, 200 000 mutations nightly under the sanitizers; the nightly suite also runs `gcc16-observe`
- [ ] OCCT differential-testing harness: build OCCT in CI, compare volumes/areas/face counts on shared cases. OCCT 8.0.1 is built and cached in CI (`tools/occt/build_occt.sh`, `.github/workflows/occt.yml`) and `stereon_occt_tests` compares the built-in curves and surfaces with OCCT's; volumes, areas and face counts follow the topology prototype
- [x] `corpus/` layout, metadata format and `corpus/LICENSES.md` ([`corpus/README.md`](corpus/README.md)); cases validated by `tools/check_corpus.py` in CI
- [ ] Regression-case workflow: failing case → minimal file + journal → never deleted
- [ ] Benchmark harness against OCCT per `docs/benchmarks.md`: pinned OCCT build in CI, Phase 0 benchmarks B0.1–B0.5 running and reported. The pinned build exists; B0.3–B0.5 run in the OCCT reference workflow with both libraries interleaved and are reported (`bench/reports/phase0-topology-reference-2026-10-20.md`); B0.2 runs there with B1 (`bench/occt/OcctBench.cpp`: direct, `VisitBuiltin` and virtual paths against `Geom_Curve::EvalD0` through `occ::handle`), its first reference report next
- [x] Seeded-fault test suite for handles and topology: stale handles, handles from another shape, off-by-one indices, wrong entity kinds, corrupted adjacency (ADR-0004 verification). `tests/topo/FaultTests.cpp`: every fault is refused through a `Result` naming the entity or by `Freeze` with its rule; the accessors' contract violations are not testable in-process

## 10. ADR prototypes (weeks 6–12)

- [ ] **ADR-0001:** union/cut/common on boxes and cylinders, including loosened-tolerance variants (1e-4 to 1e-2 mm); confirm tolerance growth is logged
- [ ] **RS-0001:** build the ADR-0001 prototype's intersections to the contract in [RS-0001](docs/research/RS-0001-tolerant-analytic-intersection.md) and run it on `corpus/analytic-ssi/`. Before that, as a side track that needs no kernel code: the literature survey, the sweep generator, the exact reference oracle and the go/no-go decision (due week 12)
- [ ] **ADR-0003:** dispatch benchmark, 10⁸ evaluations on line, circle, cubic B-spline. The `/1e8` runs and the scalar direct baseline are in `bench/GeometryBench.cpp`, and the B0.2 rows with their OCCT twins in `bench/occt/OcctBench.cpp`, run by the OCCT reference workflow; their first reference results are next
- [x] **ADR-0004:** ~10⁶-face body: adjacency speed and memory vs OCCT. Reference run of 2026-10-20 (`bench/reports/phase0-topology-reference-2026-10-20.md`): faces of every edge 20× and edges of every vertex 60× faster than OCCT; live heap per face 1255 B against OCCT's 1601 B, 78% (`docs/performance.rst`, *Topology*)
- [x] **ADR-0004:** debuggability check: seeded-fault suite caught 100% in debug builds; pretty-printers working in GDB and LLDB. The suite is in place (§9); the accessors' contract violation is tested out of process (`stereon_contract_probe shape-foreign`, `shape-range`, `shape-dead`), and the printers show the entity behind a handle in checked builds
- [x] **ADR-0004 decision point:** keep index handles throughout, or switch to the hybrid (pointers in `ShapeBuilder`, indices in `Shape`). Index handles throughout, recorded in ADR-0004 (2026-10-20)
- [ ] **ADR-0005:** 100 local edits on a 10⁶-face body (surface replacement, face deletion, face and edge splitting; no `Compact()`): memory and copy cost; page-size sweep 128/256/1 024. `bench/TopoBench.cpp` runs the sequence with replacement and deletion: 0.003% and 0.014% of the pages copied by one edit, the body and its copy after 100 edits 1.008× the body (reference run of 2026-10-20). Splitting waits for the Euler operators, and the page-size sweep is in `bench/ArenaBench.cpp`
- [ ] **ADR-0001:** merge-policy test with vertex chains at 0.8× tolerance; result deterministic and radius-bounded
- [ ] **ADR-0018:** seeded-fault corpus for the checker; zero false positives on valid bodies
- [x] **ADR-0007:** determinism check across thread counts on the prototype operations. The operations that run in parallel so far, `ProjectAll` onto B-spline curves and surfaces through `Detail::ParallelFor`, give the same bits serially and at 1, 2, 8 and 32 threads (`[determinism]`, CTest and the nightly run); booleans, tessellation and import join as they are parallelised
- [ ] **ADR-0009:** early review of the edge representation, since ADR-0004's edge record depends on it (full acceptance is due before Phase 1 topology work)
- [ ] Review ADRs 0001–0008, 0015–0018 and 0021 and move each to **Accepted** or revise. 0001, 0003, 0004, 0005, 0006 and 0016 accepted on 2026-10-20 with revisions; 0007 stays Proposed until the Boolean prototype can test it; 0018 and 0021 are open

## 11. Documentation (ongoing)

- [x] API reference generation (Doxygen or MrDocs) published to GitHub Pages. Doxygen XML through Breathe (`docs/api.rst`), built and checked by `deploy-docs.yml` (the job fails unless `api.html` documents the kernel), live at https://onurtuncer.github.io/Stereon/api.html
- [x] `docs/architecture.md`: library layering diagram and data flow. In `docs/architecture.rst`: the eleven libraries, why they stack in that order, and a Boolean from primitives to a mesh
- [ ] `docs/contributing-with-ai.md`: how agents are given ADR context and tests
- [ ] `docs/references.md`: reading list (Piegl & Tiller, Patrikalakis & Maekawa, Stroud, Golovanov, Shewchuk) with the chapters each module relies on
- [ ] `docs/debugging.md`: checked handles, `Describe()`, pretty-printer setup, journals, viewer and `stn-dump`. Started as `docs/debugging.rst`: pretty-printers, contract reports and journals; `Describe()`, the viewer and `stn-dump` follow when they exist
- [ ] Update the README roadmap status at the end of the phase

---

## Exit gate

All of the following must pass in CI before Phase 1 starts:

- [x] Predicates pass the adversarial suite with zero inconsistent results. 20 000 cases per predicate in every CI run; the 10⁶ run is `stereon_robust_tests "[exhaustive]"` (§5)
- [x] Type-erased `Curve` with `VisitBuiltin` within 10% of direct calls; pure virtual fallback within 5 ns per call of a direct call (ADR-0003, revised 2026-10-20 from 2×, a ratio that measures a line's evaluator rather than the dispatch). `VisitBuiltin` 2.0 against 2.2 ns for a line; the fallback adds 3.4 ns for a line and about 1 ns for a circle
- [x] ADR-0004 prototype: adjacency ≥ 5× faster than OCCT, memory per face < 50% of OCCT. Adjacency 20× and 60×; memory per face 608 B against OCCT's 1601 B (38%) with smaller records, per-kind wrappers and typed geometry pools (`docs/performance.rst`, *Topology*)
- [x] Release-build handles are 4 bytes; seeded-fault suite caught 100% in debug builds. `static_assert` in `src/core/Handle.cppm` and `tests/topo/RecordLayoutTests.cpp`; the faults a handle check can detect are caught by `tests/topo/FaultTests.cpp` and, at the accessors, by the contract probes
- [x] Pretty-printers pass their scripted GDB and LLDB tests. Linux Debug CI (`tests/debugger/check_printers.py`); skipped on macOS and MSVC, whose debug information the debuggers cannot read. Rerun when handles resolve to entities (§8)
- [x] Handle design decided (index throughout, or hybrid) and recorded in ADR-0004. Index throughout (2026-10-20)
- [x] Output identical for 1, 2, 8 and 32 threads. For every operation that runs in parallel today (batch projection), on every pull request and nightly across processes (§3, §10); each new parallel operation adds its case
- [x] ADR-0021 foundation types in place; no geometric operation can be added without a capability declaration
- [ ] All blocking CI jobs green on Linux, Windows and macOS
- [ ] ADRs 0001–0008, 0015–0017 accepted. All but 0007, whose speed-up and cross-platform checks wait for the Boolean prototype; 0001 accepted on 2026-10-20 with its Boolean and merge tests as conditions of that prototype
