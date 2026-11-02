# v0.1 benchmark release gates

Commit ab41d1a124182bbe874daafdaeaecef6bef3e68f, OCCT V8_0_1; verdict **FAIL**. Values are full-precision; ratios are OCCT / Stereon except the dispatch rows (Stereon over Stereon) and the counters.

| Gate | Measure | Target | Value | Threshold | Spread | Status | Note |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Run/Metadata | Commit, comparator build and host recorded | all present | — | — | — | PASS | ab41d1a124182bbe874daafdaeaecef6bef3e68f against OCCT V8_0_1 (b8f597c67) |
| Run/BenchmarkErrors | No benchmark reported an error | zero errors | — | — | — | PASS |  |
| B0.2/Dispatch/bspline | VisitBuiltin against OCCT's virtual call, cubic B-spline | faster than OCCT | 3.02123 | 1 | 0.7% | PASS | must lead outright |
| B0.2/Dispatch/circle | VisitBuiltin against OCCT's virtual call, circle | faster than OCCT | 1.19438 | 1 | 1.5% | PASS | must lead outright |
| B0.2/Dispatch/line | VisitBuiltin against OCCT's virtual call, line | faster than OCCT | 6.0046 | 1 | 0.8% | PASS | must lead outright |
| B0.3/EdgesOfVertices/1000000 | Edges of every vertex, 10⁶ faces | ≥ 5× | 52.7339 | 5 | 1.2% | PASS |  |
| B0.3/FacesOfEdges/1000000 | Faces of every edge, 10⁶ faces | ≥ 5× | 24.0169 | 5 | 26.1% | PASS |  |
| B1.1/CurveD0/cubic | Curve D0, cubic, 24 poles | ≥ parity | 3.60879 | 0.9 | 1.7% | PASS |  |
| B1.1/CurveD0/cubic-rational | Curve D0, rational cubic | ≥ parity | 2.91756 | 0.9 | 1.3% | PASS |  |
| B1.1/CurveD0/degree7 | Curve D0, degree 7 | ≥ parity | 2.21438 | 0.9 | 0.9% | PASS |  |
| B1.1/CurveD0BatchFast/cubic | Curve D0, batch of 4096, cubic, fast | ≥ 2× | 1.79477 | 2 | 1.1% | FAIL |  |
| B1.1/CurveD0BatchFast/cubic-rational | Curve D0, batch of 4096, rational cubic, fast | ≥ 2× | 1.76643 | 2 | 1.1% | FAIL |  |
| B1.1/CurveD2/cubic | Curve D0–D2, cubic | ≥ parity | 1.13297 | 0.9 | 0.6% | PASS |  |
| B1.1/CurveD2/cubic-rational | Curve D0–D2, rational cubic | ≥ parity | 1.49942 | 0.9 | 1.4% | PASS |  |
| B1.1/SurfaceD0/bicubic | Surface D0, bicubic 12 × 12 | ≥ parity | 2.49244 | 0.9 | 1.0% | PASS |  |
| B1.1/SurfaceD0/bicubic-rational | Surface D0, rational bicubic | ≥ parity | 2.69499 | 0.9 | 1.0% | PASS |  |
| B1.1/SurfaceD1/bicubic | Surface D0–D1, bicubic | ≥ parity | 1.46595 | 0.9 | 1.1% | PASS |  |
| B1.2/SurfaceGrid/256/bicubic | Grid 256 × 256, bicubic | ≥ parity | 3.17963 | 0.9 | 1.3% | PASS |  |
| B1.2/SurfaceGrid/256/bicubic-rational | Grid 256 × 256, rational bicubic | ≥ parity | 2.89314 | 0.9 | 1.6% | PASS |  |
| B1.3/CurveProjection | Curve projection, rational cubic | ≥ parity | 3.06864 | 0.9 | 7.9% | PASS |  |
| B1.3/SurfaceProjection | Surface projection, bicubic | ≥ parity | 1.10274 | 0.9 | 1.2% | PASS |  |
| B1.4/Approximation/1000 | Approximation, 1000 points, 1e-3 mm | ≥ parity at equal tolerance | 26.7197 | — | 1.5% | NOT_COMPARABLE | the libraries do not solve the same problem; a speed ratio here is no pass |
| B1.4/Interpolation/1000 | Interpolation, 1000 points, cubic | ≥ parity | 1.52298 | 0.9 | 2.1% | PASS |  |
| B1.5/BezierDecomposition | Bezier decomposition, 21 segments | ≥ parity | 1.36784 | 0.9 | 3.9% | PASS |  |
| B1.5/DegreeElevation | Degree elevation, cubic to quintic, 21 segments | ≥ parity | 0.77511 | 0.9 | 2.0% | FAIL |  |
| B1.5/KnotInsertion/100 | 100 knot insertions | ≥ parity | 2.65076 | 0.9 | 2.6% | PASS |  |
| B1.6/ArcLength/rational-cubic | Arc length of the whole curve, rational cubic | ≥ parity at equal accuracy | 1.19549 | 0.9 | 1.5% | PASS |  |
| B1.6/CurveProperties/rational-cubic | Frenet frame and curvature, rational cubic | ≥ parity | 1.40536 | 0.9 | 6.9% | PASS |  |
| B1.6/PointAtLength/64 | Parameter at a given length | ≥ parity | 1.89266 | 0.9 | 1.5% | PASS |  |
| B1.6/SurfaceProperties/bicubic | Principal curvatures and directions, bicubic | ≥ parity | 0.988261 | 0.9 | 1.2% | PASS |  |
| ADR-0004/Compact/10000 | Compact() after 50 deletions, 10⁴ faces | report | — | — | 0.5% | REPORTED |  |
| ADR-0004/Compact/1000000 | Compact() after 50 deletions, 10⁶ faces | report | — | — | 1.4% | REPORTED |  |
| ADR-0005/EditSequence/10000 | 100 edits without Compact(), 10⁴ faces | report | — | — | 1.2% | REPORTED |  |
| ADR-0005/EditSequence/1000000 | 100 edits without Compact(), 10⁶ faces | report | — | — | 0.3% | REPORTED |  |
| B0.2/Direct/bspline | Direct call, cubic B-spline, 10⁸ evaluations | report | 3.60197 | — | 0.5% | REPORTED |  |
| B0.2/Direct/circle | Direct call, circle, 10⁸ evaluations | report | 1.23216 | — | 1.5% | REPORTED |  |
| B0.2/Direct/line | Direct call, line, 10⁸ evaluations | report | 6.00253 | — | 0.3% | REPORTED |  |
| B0.2/DispatchCopy/bspline | VisitBuiltin on a local copy, cubic B-spline | report | — | — | 1.1% | REPORTED |  |
| B0.2/DispatchCopy/circle | VisitBuiltin on a local copy, circle | report | — | — | 0.8% | REPORTED |  |
| B0.2/DispatchCopy/line | VisitBuiltin on a local copy, line | report | — | — | 2.2% | REPORTED |  |
| B0.2/Virtual/bspline | Virtual call per evaluation, cubic B-spline | report | 2.9723 | — | 0.7% | REPORTED |  |
| B0.2/Virtual/circle | Virtual call per evaluation, circle | report | 1.16423 | — | 1.0% | REPORTED |  |
| B0.2/Virtual/line | Virtual call per evaluation, line | report | 4.91552 | — | 0.5% | REPORTED |  |
| B0.3/AncestorMap/10000 | OCCT's ancestor map alone, 10⁴ faces | report | — | — | 5.1% | REPORTED |  |
| B0.3/AncestorMap/1000000 | OCCT's ancestor map alone, 10⁶ faces | report | — | — | 3.4% | REPORTED |  |
| B0.3/EdgesOfVertices/10000 | Edges of every vertex, 10⁴ faces | report | 23.7035 | — | 5.3% | REPORTED |  |
| B0.3/FacesOfEdges/10000 | Faces of every edge, 10⁴ faces | report | 10.1021 | — | 3.3% | REPORTED |  |
| B0.3/SheetBuild/10000 | Build the sheet, 10⁴ faces | report | 30.8612 | — | 3.0% | REPORTED |  |
| B0.3/SheetBuild/1000000 | Build the sheet, 10⁶ faces | report | 22.8213 | — | 22.5% | REPORTED |  |
| B0.5/DeleteFace/10000 | Delete one face, 10⁴ faces | report | 10.5437 | — | 1.4% | REPORTED |  |
| B0.5/DeleteFace/1000000 | Delete one face, 10⁶ faces | report | 6.16772 | — | 1.0% | REPORTED |  |
| B0.5/ReplaceSurface/10000 | Replace one face's surface, 10⁴ faces | report | 10.6928 | — | 1.7% | REPORTED |  |
| B0.5/ReplaceSurface/1000000 | Replace one face's surface, 10⁶ faces | report | 6.29822 | — | 1.4% | REPORTED |  |
| B1.1/CurveD0Batch/cubic | Curve D0, batch of 4096, cubic, bit-identical | report | 0.515801 | — | 1.7% | REPORTED |  |
| B1.1/CurveD0Batch/cubic-rational | Curve D0, batch of 4096, rational cubic, bit-identical | report | 0.553514 | — | 1.1% | REPORTED |  |
| B0.2/DispatchOverhead/line | VisitBuiltin over the local copy, line | ≤ 1.10× | 1.00035 | 1.1 | 2.2% | PASS | VisitBuiltin / direct 1.0001 (reported) |
| B0.2/FallbackOverhead/line | Virtual call over the local copy, line | ≤ 5 ns per call | 0.149291 | 5 | 2.2% | PASS | virtual - direct 0.1491 ns (reported) |
| B0.2/DispatchOverhead/circle | VisitBuiltin over the local copy, circle | ≤ 1.10× | 0.999497 | 1.1 | 1.5% | PASS | VisitBuiltin / direct 1.0335 (reported) |
| B0.2/FallbackOverhead/circle | Virtual call over the local copy, circle | ≤ 5 ns per call | 0.274831 | 5 | 1.5% | PASS | virtual - direct 0.6337 ns (reported) |
| B0.2/DispatchOverhead/bspline | VisitBuiltin over the local copy, cubic B-spline | ≤ 1.10× | 1.00039 | 1.1 | 1.1% | PASS | VisitBuiltin / direct 1.2028 (reported) |
| B0.2/FallbackOverhead/bspline | Virtual call over the local copy, cubic B-spline | ≤ 5 ns per call | 0.380262 | 5 | 1.1% | PASS | virtual - direct 4.2219 ns (reported) |
| B0.4/HeapPerFace/1000000:heap_bytes_per_face | Live heap per face, 10⁶ faces | < 50% of OCCT | 0.379989 | 0.5 | — | PASS |  |
| B0.5/ReplaceSurface/1000000:copied_fraction | Pages copied, one surface replaced, 10⁶ faces | < 1% of body | 3.65542e-05 | 0.01 | — | PASS |  |
| B0.5/DeleteFace/1000000:copied_fraction | Pages copied, one face deleted, 10⁶ faces | < 1% of body | 9.74778e-05 | 0.01 | — | PASS |  |
| ADR-0005/EditSequence/1000000:memory_ratio | Body and its edited copy over the body, 100 edits, 10⁶ faces | < 1.2× | 1.00532 | 1.2 | — | PASS |  |
| ADR-0004/Compact/1000000:compacted_over_body | Compact() result over a fresh body, after 50 deletions, 10⁶ faces | < 1.05× | 1 | 1.05 | — | PASS |  |
| B1.3/Convergence | Projection converged on the declared corpus | ≥ 99.99% | 1 | 0.9999 | — | PASS | 20000 of 20000 converged, 0 incorrect successes, corpus b13-v1 |
