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
- [ ] One CMake target per library (`stereon::core`, `stereon::robust`, …) with private include visibility
- [x] `CMakePresets.json`: `gcc16-debug`, `gcc16-release`, `gcc16-observe`, `gcc16-asan`, `gcc16-tsan`, `gcc16-windows`, `bench`
- [ ] Contract evaluation semantics per preset: enforce (debug), observe (CI corpus), ignore/enforce-boundary (release) (ADR-0016). Debug and observe done; release is provisionally `enforce` until the per-TU split is decided
- [ ] Dependency fetching (CPM or FetchContent): stdexec, xsimd, Catch2, RapidCheck, Google Benchmark, nanobind
- [ ] Feature-test-macro wrappers: `stereon::exec` (std::execution ↔ stdexec), `stereon::simd` (std::simd ↔ xsimd)
- [ ] Include-graph check enforcing library layering (fails CI on an upward dependency)
- [ ] `.clang-format` and `.clang-tidy`, including a rule banning literal epsilon comparisons in `topo/`, `boolean/`, `mesh/` (ADR-0001, ADR-0002)
- [x] Licence-header check (REUSE or equivalent) (ADR-0008)

## 3. Continuous integration (week 2)

- [x] Linux GCC 16: debug + release, blocking
- [x] Linux GCC 16: ASan/UBSan and TSan jobs, blocking
- [x] Windows GCC 16 via MSYS2 UCRT64 (`msys2/setup-msys2`), blocking
- [x] macOS GCC 16 (Homebrew), blocking
- [ ] Clang (latest) and MSVC (latest): non-blocking, reported
- [ ] Coverage report (gcov/llvm-cov) published per PR
- [ ] Benchmark job on `main` with results stored for trend tracking
- [ ] Nightly job: full test suite + determinism check (1, 2, 8, 32 threads) (ADR-0007)

## 4. `stereon::core` (weeks 2–5)

- [ ] Fixed-size `Vec2/3`, `Point2/3`, `Mat3/4`, `Transform` with constexpr operations
- [ ] `Interval` with outward rounding
- [ ] `Box2`, `Box3` (bounding boxes) with union, intersection and interval-based tests
- [ ] Typed handles `Id<Tag>` (32-bit), one distinct type per entity kind, no implicit integer conversions (ADR-0004)
- [ ] Generational builder handles (index + generation) with free lists (ADR-0004)
- [ ] `STEREON_CHECKED_HANDLES` mode (on in debug, asan, fuzzing presets): handles carry store ID + generation; every access checks store, bounds and generation (ADR-0004 safeguards)
- [ ] `static_assert` + CI check that release-build handles are exactly 4 bytes
- [ ] Named accessors for packed orientation bits (`reversed(h)`, `base(h)`); no manual masking
- [ ] Paged arena / copy-on-write page container (`shared_ptr<const Page>`, 256 entries) (ADR-0005)
- [ ] `KernelError`, `ErrorCode`, `std::expected` aliases, `STEREON_TRY` macro (ADR-0006)
- [ ] `Context`: resolution, tolerance ceiling, scheduler, logger, journal sink (ADR-0001, ADR-0007)
- [ ] Logging and operation journal (inputs + parameters, replayable)
- [ ] Unit and property tests for all of the above

## 5. `stereon::robust` (weeks 3–6)

- [ ] Adaptive-precision `orient2d`, `orient3d`, `incircle`, `insphere` (Shewchuk-style)
- [ ] Expansion arithmetic primitives for custom predicates
- [ ] Interval-filtered predicate helpers with exact fallback
- [ ] Adversarial test suite (Shewchuk cases + 10⁶ random near-degenerate cases)
- [ ] Benchmarks: filtered vs naive floating point

## 6. Geometry scaffolding (weeks 4–8)

- [ ] Concepts `CurveGeometry`, `Curve2dGeometry`, `SurfaceGeometry` (ADR-0003)
- [ ] Type-erased `Curve`, `Curve2d`, `Surface` with 96-byte small-buffer storage and shared payloads
- [ ] Two toy types per concept: `Line` + `Circle`, `Plane` + `Cylinder`
- [ ] `visit_builtin` fast path generated with expansion statements and pack indexing (ADR-0016)
- [ ] Type-pair dispatch table for binary operations (stub intersection handlers)
- [ ] Contracts on all evaluator entry points (domain, finiteness)
- [ ] Reflection-based serialiser prototype for geometry structs (ADR-0012 groundwork)
- [ ] Property tests comparing fast path and erased path on every type

## 7. Topology prototype (weeks 6–10)

- [ ] SoA stores for vertex, edge, coedge (with radial ring), loop, face, shell, solid, following Golovanov's logical entity model (ADR-0004)
- [ ] Edge record per ADR-0009: curve kind (`Intersection`, `Exact`, `Boundary`, `Degenerate`), tolerance, exact or cached 3D curve slot; pcurves on coedges on a shared parameter
- [ ] Cache invalidation of the 3D edge curve when a pcurve changes (stub is enough for Phase 0)
- [ ] `ShapeBuilder` → `freeze()` → immutable `Shape`, with old→new remap table
- [ ] Adjacency views: `faces_of(edge)`, `edges_of(vertex)`, loop walks
- [ ] Minimal validity checks (manifoldness, loop closure, orientation), naming offending entities with `describe()`
- [ ] `describe(shape, handle)` for every entity kind, used in logs, contract messages and checker reports
- [ ] Debug builds run the validity checker after every public operation (ADR-0004 safeguards)
- [ ] Build a box (`Exact` edges) and a cylinder (`Exact` + `Boundary`/seam edges) by hand as test bodies

## 8. Debug viewer v0 (weeks 5–9)

- [ ] Qt6 or Vulkan window drawing points, polylines, boxes and triangle soups
- [ ] Load dumps written from any test (`STEREON_DUMP(obj)` macro)
- [ ] Debugger pretty-printers for GDB (Python) and LLDB: handles resolved to the entity (e.g. `Edge#4711 [Line, tol 1e-5, faces 12|37]`), `Interval`, `Box3`, `Curve`, `Surface`
- [ ] CI test for the pretty-printers using scripted GDB and LLDB sessions
- [ ] `stn-dump` prototype: JSON dump of any `Shape`, for diffing two versions (ADR-0012 groundwork)

## 9. Testing infrastructure (weeks 2–10)

- [ ] Catch2 + RapidCheck wiring; test naming and tagging conventions
- [ ] libFuzzer harness skeleton (for later STEP and primitive-placement fuzzing)
- [ ] OCCT differential-testing harness: build OCCT in CI, compare volumes/areas/face counts on shared cases
- [ ] `corpus/` layout, metadata format and `corpus/LICENSES.md`
- [ ] Regression-case workflow: failing case → minimal file + journal → never deleted
- [ ] Seeded-fault test suite for handles and topology: stale handles, handles from another shape, off-by-one indices, wrong entity kinds, corrupted adjacency (ADR-0004 verification)

## 10. ADR prototypes (weeks 6–12)

- [ ] **ADR-0001:** union/cut/common on boxes and cylinders, including loosened-tolerance variants (1e-4 to 1e-2 mm); confirm tolerance growth is logged
- [ ] **ADR-0003:** dispatch benchmark, 10⁸ evaluations on line, circle, cubic B-spline
- [ ] **ADR-0004:** ~10⁶-face body: adjacency speed and memory vs OCCT
- [ ] **ADR-0004:** debuggability check: seeded-fault suite caught 100% in debug builds; pretty-printers working in GDB and LLDB
- [ ] **ADR-0004 decision point:** keep index handles throughout, or switch to the hybrid (pointers in `ShapeBuilder`, indices in `Shape`). Decide before Phase 1 algorithms start
- [ ] **ADR-0005:** 100 local edits on a 10⁶-face body: memory and copy cost
- [ ] **ADR-0007:** determinism check across thread counts on the prototype operations
- [ ] **ADR-0009:** early review of the edge representation, since ADR-0004's edge record depends on it (full acceptance is due before Phase 1 topology work)
- [ ] Review ADRs 0001–0008, 0015, 0016 and move each to **Accepted** or revise

## 11. Documentation (ongoing)

- [ ] API reference generation (Doxygen or MrDocs) published to GitHub Pages
- [ ] `docs/architecture.md`: library layering diagram and data flow
- [ ] `docs/contributing-with-ai.md`: how agents are given ADR context and tests
- [ ] `docs/references.md`: reading list (Piegl & Tiller, Patrikalakis & Maekawa, Stroud, Golovanov, Shewchuk) with the chapters each module relies on
- [ ] `docs/debugging.md`: checked handles, `describe()`, pretty-printer setup, journals, viewer and `stn-dump`
- [ ] Update the README roadmap status at the end of the phase

---

## Exit gate

All of the following must pass in CI before Phase 1 starts:

- [ ] Predicates pass the adversarial suite with zero inconsistent results
- [ ] Type-erased `Curve` with `visit_builtin` within 10% of direct calls; pure virtual fallback within 2×
- [ ] ADR-0004 prototype: adjacency ≥ 5× faster than OCCT, memory per face < 50% of OCCT
- [ ] Release-build handles are 4 bytes; seeded-fault suite caught 100% in debug builds
- [ ] Pretty-printers pass their scripted GDB and LLDB tests
- [ ] Handle design decided (index throughout, or hybrid) and recorded in ADR-0004
- [ ] Output identical for 1, 2, 8 and 32 threads
- [ ] All blocking CI jobs green on Linux, Windows and macOS
- [ ] ADRs 0001–0008, 0015, 0016 accepted
