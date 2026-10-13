# Phase 0 — Foundations: TODO

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
- [ ] Benchmark job on `main` with results stored for trend tracking
- [ ] Nightly job: full test suite + determinism check (1, 2, 8, 32 threads) (ADR-0007)
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
- [ ] Contracts on all evaluator entry points (domain, finiteness)
- [ ] Reflection-based serialiser prototype for geometry structs (ADR-0012 groundwork)
- [x] Property tests comparing fast path and erased path on every type. Bitwise agreement for every operation of every built-in curve and surface

## 7. Topology prototype (weeks 6–10)

- [ ] SoA stores for vertex, edge, coedge (with radial ring), loop, face, shell, solid, following Golovanov's logical entity model (ADR-0004)
- [ ] Edge record per ADR-0009: curve kind (`Intersection`, `Exact`, `Trim`, `Boundary`, `Spatial`, `Degenerate`), tolerance, stored primary coedge, own or cached 3D curve slot with certificate; pcurves on coedges on a shared parameter
- [ ] Cache invalidation of the 3D edge curve when a pcurve changes (stub is enough for Phase 0)
- [ ] `ShapeBuilder` → `Freeze()` → immutable `Shape` with stable indices and alive bitsets; separate `Compact()` returning the old→new remap table (ADR-0004)
- [ ] Adjacency views: `FacesOf(edge)`, `EdgesOf(vertex)`, loop walks
- [ ] `Structural`-level validity checker with rule IDs and a `ValidityReport`, naming offending entities with `Describe()` (ADR-0018); `Geometric` level once edge certificates exist
- [ ] `OperationResult` (shape, history, tolerance changes, warnings) as the return type of shape-producing operations (ADR-0006)
- [ ] `Describe(shape, handle)` for every entity kind, used in logs, contract messages and checker reports
- [ ] Debug builds run the validity checker after every public operation (ADR-0004 safeguards)
- [ ] Build a box (`Exact` edges) and a cylinder (`Exact` + `Boundary`/seam edges) by hand as test bodies

## 8. Debug viewer v0 (weeks 5–9)

- [ ] Qt6 or Vulkan window drawing points, polylines, boxes and triangle soups
- [ ] Load dumps written from any test (`STEREON_DUMP(obj)` macro)
- [ ] Debugger pretty-printers for GDB (Python) and LLDB: handles resolved to the entity (e.g. `Edge#4711 [Line, tol 1e-5, faces 12|37]`), `Interval`, `Box3`, `Curve`, `Surface`
- [ ] CI test for the pretty-printers using scripted GDB and LLDB sessions
- [ ] `stn-dump` prototype: JSON dump of any `Shape`, for diffing two versions (ADR-0012 groundwork)

## 9. Testing infrastructure (weeks 2–10)

- [x] Catch2 + RapidCheck wiring; test naming and tagging conventions. One executable per library, cases tagged `[<library>]` plus topic tags and `[property]` for `rc::prop` tests ([`tests/CMakeLists.txt`](tests/CMakeLists.txt))
- [ ] libFuzzer harness skeleton (for later STEP and primitive-placement fuzzing)
- [ ] OCCT differential-testing harness: build OCCT in CI, compare volumes/areas/face counts on shared cases
- [x] `corpus/` layout, metadata format and `corpus/LICENSES.md` ([`corpus/README.md`](corpus/README.md)); cases validated by `tools/check_corpus.py` in CI
- [ ] Regression-case workflow: failing case → minimal file + journal → never deleted
- [ ] Benchmark harness against OCCT per `docs/benchmarks.md`: pinned OCCT build in CI, Phase 0 benchmarks B0.1–B0.5 running and reported
- [ ] Seeded-fault test suite for handles and topology: stale handles, handles from another shape, off-by-one indices, wrong entity kinds, corrupted adjacency (ADR-0004 verification)

## 10. ADR prototypes (weeks 6–12)

- [ ] **ADR-0001:** union/cut/common on boxes and cylinders, including loosened-tolerance variants (1e-4 to 1e-2 mm); confirm tolerance growth is logged
- [ ] **RS-0001:** build the ADR-0001 prototype's intersections to the contract in [RS-0001](docs/research/RS-0001-tolerant-analytic-intersection.md) and run it on `corpus/analytic-ssi/`. Before that, as a side track that needs no kernel code: the literature survey, the sweep generator, the exact reference oracle and the go/no-go decision (due week 12)
- [ ] **ADR-0003:** dispatch benchmark, 10⁸ evaluations on line, circle, cubic B-spline
- [ ] **ADR-0004:** ~10⁶-face body: adjacency speed and memory vs OCCT
- [ ] **ADR-0004:** debuggability check: seeded-fault suite caught 100% in debug builds; pretty-printers working in GDB and LLDB
- [ ] **ADR-0004 decision point:** keep index handles throughout, or switch to the hybrid (pointers in `ShapeBuilder`, indices in `Shape`). Decide before Phase 1 algorithms start
- [ ] **ADR-0005:** 100 local edits on a 10⁶-face body (surface replacement, face deletion, face and edge splitting; no `Compact()`): memory and copy cost; page-size sweep 128/256/1 024
- [ ] **ADR-0001:** merge-policy test with vertex chains at 0.8× tolerance; result deterministic and radius-bounded
- [ ] **ADR-0018:** seeded-fault corpus for the checker; zero false positives on valid bodies
- [ ] **ADR-0007:** determinism check across thread counts on the prototype operations
- [ ] **ADR-0009:** early review of the edge representation, since ADR-0004's edge record depends on it (full acceptance is due before Phase 1 topology work)
- [ ] Review ADRs 0001–0008, 0015–0018 and move each to **Accepted** or revise

## 11. Documentation (ongoing)

- [ ] API reference generation (Doxygen or MrDocs) published to GitHub Pages
- [ ] `docs/architecture.md`: library layering diagram and data flow
- [ ] `docs/contributing-with-ai.md`: how agents are given ADR context and tests
- [ ] `docs/references.md`: reading list (Piegl & Tiller, Patrikalakis & Maekawa, Stroud, Golovanov, Shewchuk) with the chapters each module relies on
- [ ] `docs/debugging.md`: checked handles, `Describe()`, pretty-printer setup, journals, viewer and `stn-dump`
- [ ] Update the README roadmap status at the end of the phase

---

## Exit gate

All of the following must pass in CI before Phase 1 starts:

- [ ] Predicates pass the adversarial suite with zero inconsistent results
- [ ] Type-erased `Curve` with `VisitBuiltin` within 10% of direct calls; pure virtual fallback within 2×
- [ ] ADR-0004 prototype: adjacency ≥ 5× faster than OCCT, memory per face < 50% of OCCT
- [ ] Release-build handles are 4 bytes; seeded-fault suite caught 100% in debug builds
- [ ] Pretty-printers pass their scripted GDB and LLDB tests
- [ ] Handle design decided (index throughout, or hybrid) and recorded in ADR-0004
- [ ] Output identical for 1, 2, 8 and 32 threads
- [ ] All blocking CI jobs green on Linux, Windows and macOS
- [ ] ADRs 0001–0008, 0015–0017 accepted
