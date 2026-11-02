# v0.1 benchmark release gates

Commit 1449128db71a50606f77bd290bd5fc3b00acf798, OCCT V8_0_1; verdict **FAIL**. Values are full-precision; ratios are OCCT / Stereon except the dispatch rows (Stereon over Stereon) and the counters.

| Gate | Measure | Target | Value | Threshold | Spread | Status | Note |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Run/Metadata | Commit, comparator build and host recorded | all present | — | — | — | PASS | 1449128db71a50606f77bd290bd5fc3b00acf798 against OCCT V8_0_1 (b8f597c67) |
| Run/BenchmarkErrors | No benchmark reported an error | zero errors | — | — | — | PASS |  |
| B0.2/Dispatch/bspline | VisitBuiltin against OCCT's virtual call, cubic B-spline | faster than OCCT | 2.84955 | 1 | 17.9% | PASS | must lead outright |
| B0.2/Dispatch/circle | VisitBuiltin against OCCT's virtual call, circle | faster than OCCT | 1.22673 | 1 | 10.0% | PASS | must lead outright |
| B0.2/Dispatch/line | VisitBuiltin against OCCT's virtual call, line | faster than OCCT | 4.99918 | 1 | 4.4% | PASS | must lead outright |
| B0.3/EdgesOfVertices/1000000 | Edges of every vertex, 10⁶ faces | ≥ 5× | 69.2138 | 5 | 21.7% | PASS |  |
| B0.3/FacesOfEdges/1000000 | Faces of every edge, 10⁶ faces | ≥ 5× | 25.9218 | 5 | 2.9% | PASS |  |
| B1.1/CurveD0/cubic | Curve D0, cubic, 24 poles | ≥ parity | 3.8584 | 0.9 | 11.4% | PASS |  |
| B1.1/CurveD0/cubic-rational | Curve D0, rational cubic | ≥ parity | 3.00243 | 0.9 | 3.4% | PASS |  |
| B1.1/CurveD0/degree7 | Curve D0, degree 7 | ≥ parity | 2.22978 | 0.9 | 7.6% | PASS |  |
| B1.1/CurveD0BatchFast/cubic | Curve D0, batch of 4096, cubic, fast | ≥ 2× | 2.0032 | 2 | 6.4% | PASS |  |
| B1.1/CurveD0BatchFast/cubic-rational | Curve D0, batch of 4096, rational cubic, fast | ≥ 2× | 1.99379 | 2 | 0.9% | FAIL |  |
| B1.1/CurveD2/cubic | Curve D0–D2, cubic | ≥ parity | 0.960488 | 0.9 | 6.6% | PASS |  |
| B1.1/CurveD2/cubic-rational | Curve D0–D2, rational cubic | ≥ parity | 1.30924 | 0.9 | 8.9% | PASS |  |
| B1.1/SurfaceD0/bicubic | Surface D0, bicubic 12 × 12 | ≥ parity | 2.91348 | 0.9 | 4.3% | PASS |  |
| B1.1/SurfaceD0/bicubic-rational | Surface D0, rational bicubic | ≥ parity | 2.8297 | 0.9 | 2.4% | PASS |  |
| B1.1/SurfaceD1/bicubic | Surface D0–D1, bicubic | ≥ parity | 1.21458 | 0.9 | 1.5% | PASS |  |
| B1.2/SurfaceGrid/256/bicubic | Grid 256 × 256, bicubic | ≥ parity | 3.5552 | 0.9 | 3.1% | PASS |  |
| B1.2/SurfaceGrid/256/bicubic-rational | Grid 256 × 256, rational bicubic | ≥ parity | 3.48862 | 0.9 | 0.3% | PASS |  |
| B1.3/CurveProjection | Curve projection, rational cubic | ≥ parity | 3.39813 | 0.9 | 1.1% | PASS |  |
| B1.3/SurfaceProjection | Surface projection, bicubic | ≥ parity | 1.07823 | 0.9 | 1.5% | PASS |  |
| B1.4/Approximation/1000 | Approximation, 1000 points, 1e-3 mm | ≥ parity at equal tolerance | 23.2658 | — | 1.5% | NOT_COMPARABLE | the libraries do not solve the same problem; a speed ratio here is no pass |
| B1.4/Interpolation/1000 | Interpolation, 1000 points, cubic | ≥ parity | 1.6773 | 0.9 | 1.3% | PASS |  |
| B1.5/BezierDecomposition | Bezier decomposition, 21 segments | ≥ parity | 1.36028 | 0.9 | 6.6% | PASS |  |
| B1.5/DegreeElevation | Degree elevation, cubic to quintic, 21 segments | ≥ parity | 1.05469 | 0.9 | 6.5% | PASS |  |
| B1.5/KnotInsertion/100 | 100 knot insertions | ≥ parity | 2.39299 | 0.9 | 1.0% | PASS |  |
| B1.6/ArcLength/rational-cubic | Arc length of the whole curve, rational cubic | ≥ parity at equal accuracy | 0.945486 | 0.9 | 5.4% | PASS |  |
| B1.6/CurveProperties/rational-cubic | Frenet frame and curvature, rational cubic | ≥ parity | 1.30321 | 0.9 | 5.7% | PASS |  |
| B1.6/PointAtLength/64 | Parameter at a given length | ≥ parity | 1.51488 | 0.9 | 12.1% | PASS |  |
| B1.6/SurfaceProperties/bicubic | Principal curvatures and directions, bicubic | ≥ parity | 0.898816 | 0.9 | 3.9% | FAIL |  |
| ADR-0004/Compact/10000 | Compact() after 50 deletions, 10⁴ faces | report | — | — | 3.1% | REPORTED |  |
| ADR-0004/Compact/1000000 | Compact() after 50 deletions, 10⁶ faces | report | — | — | 2.3% | REPORTED |  |
| ADR-0005/EditSequence/10000 | 100 edits without Compact(), 10⁴ faces | report | — | — | 4.2% | REPORTED |  |
| ADR-0005/EditSequence/1000000 | 100 edits without Compact(), 10⁶ faces | report | — | — | 0.8% | REPORTED |  |
| B0.2/Direct/bspline | Direct call, cubic B-spline, 10⁸ evaluations | report | 3.5644 | — | 4.4% | REPORTED |  |
| B0.2/Direct/circle | Direct call, circle, 10⁸ evaluations | report | 1.28195 | — | 10.1% | REPORTED |  |
| B0.2/Direct/line | Direct call, line, 10⁸ evaluations | report | 5.00353 | — | 33.4% | REPORTED |  |
| B0.2/DispatchCopy/bspline | VisitBuiltin on a local copy, cubic B-spline | report | — | — | 1.4% | REPORTED |  |
| B0.2/DispatchCopy/circle | VisitBuiltin on a local copy, circle | report | — | — | 2.0% | REPORTED |  |
| B0.2/DispatchCopy/line | VisitBuiltin on a local copy, line | report | — | — | 0.1% | REPORTED |  |
| B0.2/Virtual/bspline | Virtual call per evaluation, cubic B-spline | report | 3.27181 | — | 14.6% | REPORTED |  |
| B0.2/Virtual/circle | Virtual call per evaluation, circle | report | 1.23671 | — | 11.2% | REPORTED |  |
| B0.2/Virtual/line | Virtual call per evaluation, line | report | 4.99536 | — | 19.1% | REPORTED |  |
| B0.3/AncestorMap/10000 | OCCT's ancestor map alone, 10⁴ faces | report | — | — | 7.3% | REPORTED |  |
| B0.3/AncestorMap/1000000 | OCCT's ancestor map alone, 10⁶ faces | report | — | — | 3.1% | REPORTED |  |
| B0.3/EdgesOfVertices/10000 | Edges of every vertex, 10⁴ faces | report | 30.1345 | — | 7.6% | REPORTED |  |
| B0.3/FacesOfEdges/10000 | Faces of every edge, 10⁴ faces | report | 8.99972 | — | 5.1% | REPORTED |  |
| B0.3/SheetBuild/10000 | Build the sheet, 10⁴ faces | report | 23.4433 | — | 2.7% | REPORTED |  |
| B0.3/SheetBuild/1000000 | Build the sheet, 10⁶ faces | report | 17.9514 | — | 6.4% | REPORTED |  |
| B0.5/DeleteFace/10000 | Delete one face, 10⁴ faces | report | 5.71379 | — | 13.9% | REPORTED |  |
| B0.5/DeleteFace/1000000 | Delete one face, 10⁶ faces | report | 3.66332 | — | 1.6% | REPORTED |  |
| B0.5/ReplaceSurface/10000 | Replace one face's surface, 10⁴ faces | report | 5.86356 | — | 11.0% | REPORTED |  |
| B0.5/ReplaceSurface/1000000 | Replace one face's surface, 10⁶ faces | report | 3.77495 | — | 8.7% | REPORTED |  |
| B1.1/CurveD0Batch/cubic | Curve D0, batch of 4096, cubic, bit-identical | report | 0.358318 | — | 4.7% | REPORTED |  |
| B1.1/CurveD0Batch/cubic-rational | Curve D0, batch of 4096, rational cubic, bit-identical | report | 0.421434 | — | 0.7% | REPORTED |  |
| B0.2/DispatchOverhead/line | VisitBuiltin over the local copy, line | ≤ 1.10× | 1.00014 | 1.1 | 19.1% | INCONCLUSIVE | VisitBuiltin / direct 1.0008 (reported) |
| B0.2/FallbackOverhead/line | Virtual call over the local copy, line | ≤ 5 ns per call | 0.000871531 | 5 | 19.1% | PASS | virtual - direct 0.0014 ns (reported) |
| B0.2/DispatchOverhead/circle | VisitBuiltin over the local copy, circle | ≤ 1.10× | 0.999765 | 1.1 | 2.0% | PASS | VisitBuiltin / direct 1.0551 (reported) |
| B0.2/FallbackOverhead/circle | Virtual call over the local copy, circle | ≤ 5 ns per call | -0.0855293 | 5 | 2.0% | PASS | virtual - direct 0.5558 ns (reported) |
| B0.2/DispatchOverhead/bspline | VisitBuiltin over the local copy, cubic B-spline | ≤ 1.10× | 0.998489 | 1.1 | 1.5% | PASS | VisitBuiltin / direct 1.2021 (reported) |
| B0.2/FallbackOverhead/bspline | Virtual call over the local copy, cubic B-spline | ≤ 5 ns per call | 0.415481 | 5 | 1.5% | PASS | virtual - direct 4.1576 ns (reported) |
| B0.4/HeapPerFace/1000000:heap_bytes_per_face | Live heap per face, 10⁶ faces | < 50% of OCCT | 0.379989 | 0.5 | — | PASS |  |
| B0.5/ReplaceSurface/1000000:copied_fraction | Pages copied, one surface replaced, 10⁶ faces | < 1% of body | 3.65542e-05 | 0.01 | — | PASS |  |
| B0.5/DeleteFace/1000000:copied_fraction | Pages copied, one face deleted, 10⁶ faces | < 1% of body | 9.74778e-05 | 0.01 | — | PASS |  |
| ADR-0005/EditSequence/1000000:memory_ratio | Body and its edited copy over the body, 100 edits, 10⁶ faces | < 1.2× | 1.00532 | 1.2 | — | PASS |  |
| ADR-0004/Compact/1000000:compacted_over_body | Compact() result over a fresh body, after 50 deletions, 10⁶ faces | < 1.05× | 1 | 1.05 | — | PASS |  |
| B1.3/Convergence | Projection converged on the declared corpus | ≥ 99.99% | 1 | 0.9999 | — | PASS | 20000 of 20000 converged, 0 incorrect successes, corpus b13-v1 |
