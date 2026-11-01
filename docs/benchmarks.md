# Benchmarks against OpenCASCADE

Stereon is measured against OpenCASCADE (OCCT) on **speed, memory and robustness**, phase by phase. A speed result only counts when both kernels produce a valid result for the same case; failures are reported separately as a robustness score.

## Methodology

- **Same machine, pinned versions.** OCCT built from a pinned tag in Release mode with the same compiler (GCC 16) and flags (`-O2`, LTO off for both, same `-march`). Versions are recorded in every report.
  The pin lives in `tools/occt/build_occt.sh` (OCCT `V8_0_1`, commit `b8f597c`, modelling modules only), which CI runs and caches; each install records its tag, commit and compiler in `stereon-occt-build.txt`. The OCCT build uses CMake's Release flags; matching them to Stereon's benchmark flags is part of the benchmark harness.
- **Fair parallelism.** Each benchmark runs single-threaded and at 8 threads. OCCT's parallel switches are enabled where they exist (`BOPAlgo_Options::SetRunParallel`, `BRepMesh_IncrementalMesh` parallel flag, `OSD_ThreadPool`).
- **Repetitions and spread.** The reference workflow (`.github/workflows/occt-reference.yml`) runs each curve and surface benchmark 15 times and each topology benchmark 5 times, both libraries interleaved at random on the same runner, and reports the median with the spread (the largest deviation of any repetition from its median); a row whose spread is large is re-run, not quoted. Cold (first-run) time is reported separately for I/O benchmarks. *Revised 2026-10-31 and approved by the owner the same day:* this line first asked for the median and p95 of ≥ 20 warm runs, which no reference run produced; the committed reports record 5 and 15 repetitions with the spread. The revision documents the protocol in use, now the protocol of record. On a shared hosted runner the host's variation between runs, up to about 1.5× in absolute time on both sides, bounds the precision, not the repetition count, so ratios are read within one run and within 10 % of 1 as parity.
- **Correctness gate.** Every result is checked: Stereon with its validity checker, OCCT with `BRepCheck_Analyzer`. Volumes/areas compared with `BRepGProp` and analytic values. Mismatches > 1e-6 relative are flagged, not averaged in.
- **Robustness score.** For each corpus: % valid results, % structured errors (Stereon `KernelError`, OCCT `IsDone() == false` or exception), % *silently invalid* results. The last must be 0% for Stereon.
- **Memory.** Peak RSS and retained bytes per entity (heap profiler), measured after the operation.
- **Harness.** Google Benchmark for micro-benchmarks; a corpus runner for operation-level benchmarks; results stored per commit and plotted as trends. OCCT↔Stereon geometry is exchanged through the native-format bridge (see ADR-0004 prior art: `BRepGraph` tables can feed it).

## Datasets

| Corpus | Used for |
| --- | --- |
| Synthetic generators (random primitives, patterned arrays, near-tangent and coincident placements) | Booleans, classification, stress at scale |
| SSI corpus (≥ 500 surface pairs, ADR-0010) | Intersection speed and completeness |
| ABC dataset sample, Fusion 360 Gallery | Real parts: booleans, tessellation, healing |
| NIST CAD test models | STEP import/export |
| Large bodies (10⁵–10⁶ faces) | Adjacency, memory, snapshots, tessellation scaling |

## Benchmark list

### Phase 0 — foundations

| ID | Benchmark | OCCT counterpart | Metric | Target |
| --- | --- | --- | --- | --- |
| B0.1 | Robust `Orient3D` on random and near-degenerate points | Plain floating-point determinant (OCCT has no exact predicate) | ns/call; inconsistency count | ≤ 1.5× naive cost, 0 inconsistencies |
| B0.2 | Curve dispatch: 10⁸ evaluations, line/circle/cubic B-spline | `Geom_Curve::EvalD0` through `occ::handle` | ns/eval | `VisitBuiltin` ≤ 1.1× the same loop through `VisitBuiltin` on a local copy of the geometry (`B0.2/DispatchCopy`; ADR-0003, revised 2026-11-01), and faster than OCCT's virtual call; the ratio to the direct call on the concrete type is reported, not gated; the virtual fallback adds ≤ 5 ns per call over the direct call |
| B0.3 | Adjacency on a 10⁶-face body: faces of each edge, edges of each vertex | `TopExp::MapShapesAndAncestors`; also `BRepGraph` relation queries | total time; memory | ≥ 5× faster than `TopExp`; report vs `BRepGraph` |
| B0.4 | Memory per entity | Same body as `TopoDS` and as `BRepGraph` | bytes/face, bytes/edge | < 50% of `TopoDS` |
| B0.5 | Copy / snapshot after editing one face of a 10⁶-face body | `BRepBuilderAPI_Copy`; `BRepGraph_Copy` | time; bytes copied | < 1% of body copied |

### Phase 1 — curves and surfaces

| ID | Benchmark | OCCT counterpart | Metric | Target |
| --- | --- | --- | --- | --- |
| B1.1 | NURBS curve/surface evaluation (D0–D2), degrees 2–7 | `Geom_BSplineSurface::D0/D1/D2` | ns/eval | ≥ parity; fast SIMD batch (`EvalAllFast`) ≥ 2× `GeomGridEval_BSplineCurve`; the bit-identical batch (`EvalAll`) is reported without a target |
| B1.2 | Grid evaluation (e.g. 256×256) | `GeomGridEval_BSplineSurface` | points/s | ≥ parity |
| B1.3 | Point projection / inversion on NURBS surfaces | `GeomAPI_ProjectPointOnSurf`, `Extrema_ExtPS` | time; % converged; max error | ≥ parity, 99.99% convergence |
| B1.4 | Curve and surface fitting / interpolation | `GeomAPI_PointsToBSpline(Surface)`, `GeomAPI_Interpolate` | time; max deviation | ≥ parity at equal tolerance |
| B1.5 | Knot insertion, degree elevation, Bézier decomposition | `Geom_BSplineCurve::InsertKnot`, `IncreaseDegree`, `GeomConvert` | time | ≥ parity |
| B1.6 | Local properties (Frenet frame, principal curvatures), arc length, point at a length | `GeomLProp_CLProps`, `GeomLProp_SLProps`, `GCPnts_AbscissaPoint` | time; length agreement | ≥ parity at equal accuracy |

### Phase 2 — intersections

| ID | Benchmark | OCCT counterpart | Metric | Target |
| --- | --- | --- | --- | --- |
| B2.1 | Curve–curve, curve–surface intersection | `Geom2dAPI_InterCurveCurve`, `GeomAPI_IntCS` | time; branches found | parity speed, no missed roots |
| B2.2 | Surface–surface intersection on the SSI corpus | `GeomAPI_IntSS` / `IntPatch` | time; missed / spurious branches; max deviation | 0 missed, 0 spurious; faster on median |
| B2.3 | Analytic SSI (plane, quadric pairs) | `IntAna_QuadQuadGeo` | time; exactness | exact results, ≥ parity |

### Phase 3 — topology and construction

| ID | Benchmark | OCCT counterpart | Metric | Target |
| --- | --- | --- | --- | --- |
| B3.1 | Primitives, extrude, revolve, sweep, loft | `BRepPrimAPI_*`, `BRepOffsetAPI_MakePipeShell`, `BRepOffsetAPI_ThruSections` | time; validity | ≥ parity, 100% valid |
| B3.2 | Validity check of large bodies | `BRepCheck_Analyzer` | time | ≥ 3× faster |
| B3.3 | Mass properties (volume, area, centroid, inertia) | `BRepGProp::VolumeProperties` | time; error vs analytic | ≥ parity, 1e-10 relative |
| B3.4 | Sewing a face soup | `BRepBuilderAPI_Sewing` | time; % closed shells | ≥ parity |
| B3.5 | Native save / load of 10⁶-face body | `BinTools` / `BRepTools::Write` | MB/s; file size | load < 1 s, ≥ 3× faster |

### Phase 4 — booleans (the headline benchmark)

| ID | Benchmark | OCCT counterpart | Metric | Target |
| --- | --- | --- | --- | --- |
| B4.1 | Union/cut/common on the 5 000-case corpus | `BRepAlgoAPI_Fuse/Cut/Common` | robustness score; median and p95 time | ≥ 99.5% valid, 0% silently invalid, faster median |
| B4.2 | Degenerate cases: coincident faces, tangent contact, shared edges | same | robustness score | ≥ 99% valid |
| B4.3 | Patterned arrays (e.g. plate minus 10⁴ holes) | same, with `SetRunParallel(true)` | time vs count; 1 and 8 threads | ≥ 5× speed-up at 8 threads |
| B4.4 | Point-in-solid classification | `BRepClass3d_SolidClassifier` | queries/s; errors | ≥ parity, 0 misclassifications |
| B4.5 | Tolerance growth after booleans | same | max and mean edge/vertex tolerance after N chained booleans | no unlogged growth; lower than OCCT |
| B4.6 | Determinism | same | output identical across thread counts | identical (OCCT reported for information) |

### Phase 5 — tessellation and exchange

| ID | Benchmark | OCCT counterpart | Metric | Target |
| --- | --- | --- | --- | --- |
| B5.1 | Display tessellation at fixed chordal/angular tolerance | `BRepMesh_IncrementalMesh` | triangles/s; triangle count; max chordal error | ≥ parity speed, ≤ triangle count at equal error |
| B5.2 | Analysis-grade watertight mesh | `BRepMesh` + check | % watertight; quality (min angle) | 100% watertight |
| B5.3 | STEP AP242 import | `STEPControl_Reader` (+ `ShapeFix_Shape`) | MB/s; % valid solids | ≥ parity speed, ≥ 95% valid |
| B5.4 | STEP export | `STEPControl_Writer` | MB/s; round-trip fidelity | ≥ parity |
| B5.5 | Healing | `ShapeFix_Shape`, `ShapeUpgrade` | % repaired; tolerance after repair | ≥ parity |
| B5.6 | Minimum distance between shapes | `BRepExtrema_DistShapeShape` | time; error | ≥ parity |
| B5.7 | B-rep to tri-dexel conversion (ADR-0044) on identical grids | `IntCurvesFace_ShapeIntersector` per ray; `BRepMesh` + ray casting for information | rays/s at 1 and 8 threads; peak memory; endpoint error; unresolved and silently wrong rays | ≥ 3× at 8 threads at equal endpoint accuracy, 0 silently wrong |

### Phase 6 — blends and advanced modelling

| ID | Benchmark | OCCT counterpart | Metric | Target |
| --- | --- | --- | --- | --- |
| B6.1 | Constant-radius fillets on edge chains | `BRepFilletAPI_MakeFillet` | robustness score; time | ≥ 95% valid, ≥ OCCT success rate |
| B6.2 | Chamfers | `BRepFilletAPI_MakeChamfer` | robustness score; time | ≥ parity |
| B6.3 | Shelling / thickening | `BRepOffsetAPI_MakeThickSolid` | robustness score; time | ≥ OCCT success rate |
| B6.4 | Shape sensitivities | none in OCCT (finite differences used as baseline) | time vs finite differences; accuracy | ≥ 5× faster than FD, matches FD |

## Reporting

- One report per release gate, stored in `bench/reports/`, with OCCT version, hardware, and the full table above.
- The report of a release gate is a **reference run**: the `OCCT reference run` workflow (`.github/workflows/occt-reference.yml`) on the reference image, `gcc:16.2` on a GitHub-hosted `ubuntu-24.04` runner, against the cached OCCT pin. It uploads the raw JSON, the machine record and the report that `tools/bench_report.py` renders; the release commits that report. Runs on a developer's machine are development runs and say so.
- The trend between releases is kept by the `Benchmarks` workflow (`.github/workflows/bench.yml`): on every push to main, `stereon_bench` with five repetitions, recorded by `tools/bench_trend.py` on the `bench-data` branch, compared with the run before and plotted on the trend page published with the manual (`benchmarks/trend/`). Stereon only, and never a gate: the runner differs between runs.
- Release notes quote only benchmarks whose correctness gate passed on both kernels.
- A regression of > 10% in any Stereon benchmark between commits blocks the merge unless explained in the PR.

### M5 splitting edits

`B0.5/SplitEdge/Stereon` and `B0.5/SplitFace/Stereon` use the existing
10,000- and 1,000,000-face sheet fixtures and report elapsed time, copied
pages, copied fraction and heap delta through the same `OneEdit` harness.
The edge case splits a shared interior edge and all its radial uses; the face
case inserts a diagonal. Run these in the native benchmark configuration
before accepting ADR-0005's copying targets. No measurements are recorded yet.

`ADR-0005/EditSequence/Stereon` now includes 25 edits of each of the four
kinds, on distinct original faces, and reports the maximum per-edit copied
fraction as well as the combined snapshot memory ratio. The earlier two-kind
measurements do not qualify this expanded sequence.
