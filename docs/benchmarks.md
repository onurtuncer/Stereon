# Benchmarks against OpenCASCADE

Stereon is measured against OpenCASCADE (OCCT) on **speed, memory and robustness**, phase by phase. A speed result only counts when both kernels produce a valid result for the same case; failures are reported separately as a robustness score.

## Methodology

- **Same machine, pinned versions.** OCCT built from a pinned tag in Release mode with the same compiler (GCC 16) and flags (`-O2`, LTO off for both, same `-march`). Versions are recorded in every report.
- **Fair parallelism.** Each benchmark runs single-threaded and at 8 threads. OCCT's parallel switches are enabled where they exist (`BOPAlgo_Options::SetRunParallel`, `BRepMesh_IncrementalMesh` parallel flag, `OSD_ThreadPool`).
- **Warm and cold.** Report the median and p95 of ≥ 20 warm runs; cold (first-run) time is reported separately for I/O benchmarks.
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
| B0.2 | Curve dispatch: 10⁸ evaluations, line/circle/cubic B-spline | `Geom_Curve::EvalD0` through `occ::handle` | ns/eval | `VisitBuiltin` ≤ 1.1× direct call and faster than OCCT virtual call |
| B0.3 | Adjacency on a 10⁶-face body: faces of each edge, edges of each vertex | `TopExp::MapShapesAndAncestors`; also `BRepGraph` relation queries | total time; memory | ≥ 5× faster than `TopExp`; report vs `BRepGraph` |
| B0.4 | Memory per entity | Same body as `TopoDS` and as `BRepGraph` | bytes/face, bytes/edge | < 50% of `TopoDS` |
| B0.5 | Copy / snapshot after editing one face of a 10⁶-face body | `BRepBuilderAPI_Copy`; `BRepGraph_Copy` | time; bytes copied | < 1% of body copied |

### Phase 1 — curves and surfaces

| ID | Benchmark | OCCT counterpart | Metric | Target |
| --- | --- | --- | --- | --- |
| B1.1 | NURBS curve/surface evaluation (D0–D2), degrees 2–7 | `Geom_BSplineSurface::D0/D1/D2` | ns/eval | ≥ parity; SIMD batch ≥ 2× |
| B1.2 | Grid evaluation (e.g. 256×256) | `GeomGridEval_BSplineSurface` | points/s | ≥ parity |
| B1.3 | Point projection / inversion on NURBS surfaces | `GeomAPI_ProjectPointOnSurf`, `Extrema_ExtPS` | time; % converged; max error | ≥ parity, 99.99% convergence |
| B1.4 | Curve and surface fitting / interpolation | `GeomAPI_PointsToBSpline(Surface)`, `GeomAPI_Interpolate` | time; max deviation | ≥ parity at equal tolerance |
| B1.5 | Knot insertion, degree elevation, Bézier decomposition | `Geom_BSplineCurve::InsertKnot`, `IncreaseDegree`, `GeomConvert` | time | ≥ parity |

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

### Phase 6 — blends and advanced modelling

| ID | Benchmark | OCCT counterpart | Metric | Target |
| --- | --- | --- | --- | --- |
| B6.1 | Constant-radius fillets on edge chains | `BRepFilletAPI_MakeFillet` | robustness score; time | ≥ 95% valid, ≥ OCCT success rate |
| B6.2 | Chamfers | `BRepFilletAPI_MakeChamfer` | robustness score; time | ≥ parity |
| B6.3 | Shelling / thickening | `BRepOffsetAPI_MakeThickSolid` | robustness score; time | ≥ OCCT success rate |
| B6.4 | Shape sensitivities | none in OCCT (finite differences used as baseline) | time vs finite differences; accuracy | ≥ 5× faster than FD, matches FD |

## Reporting

- One report per release gate, stored in `bench/reports/`, with OCCT version, hardware, and the full table above.
- Release notes quote only benchmarks whose correctness gate passed on both kernels.
- A regression of > 10% in any Stereon benchmark between commits blocks the merge unless explained in the PR.
