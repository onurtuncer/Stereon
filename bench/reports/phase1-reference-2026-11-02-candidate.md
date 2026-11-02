# Benchmark reference run — dispatch (B0.2), topology (B0.3–B0.5), curves and surfaces (B1.1–B1.6)

## Setup

| | |
| --- | --- |
| Date | 2026-11-02 |
| Commit | 1449128db71a50606f77bd290bd5fc3b00acf798 |
| Ref | main |
| Workflow run | https://github.com/onurtuncer/Stereon/actions/runs/41460843836 |
| Runner | GitHub Actions ubuntu-24.04 (Linux, X64), container gcc:16.2 |
| CPU | AMD EPYC 9V74 80-Core Processor |
| CPUs | 4 |
| Memory | 15.6 GB |
| Kernel | Linux 7.0.0-1017-azure |
| Compiler | g++ (GCC) 16.2.0 |
| OCCT tag | V8_0_1 |
| OCCT commit | b8f597c677811d1f9f4d8a97f5ae2825c0353a42 |
| OCCT compiler | c++ (GCC) 16.2.0 |
| OCCT build type | Release |
| Repetitions | 5, 15 |
| memory_measure | live heap, glibc mallinfo2 (bench/Memory.hpp) |
| topology | threads per core 2,cores per socket 2,sockets 1 |

## Results

Median CPU time over the repetitions, both libraries interleaved: per evaluated parameter, projected point, property or length query; per fit (B1.4), per 100 insertions and per decomposition (B1.5); per build, per query over the whole sheet or per edit (B0.3, B0.5, ADR-0004, ADR-0005). Ratio = OCCT / Stereon, so > 1 means Stereon is faster; within 10 % of 1 is parity, and a row with a ratio target of its own (B0.3: ≥ 5×) is met when it reaches it. A target of "report" is not a gate. Spread is the largest relative deviation of any repetition from its median, over both libraries.

| ID | Benchmark | Stereon | OCCT | Ratio | Spread | Target | Met |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ADR-0004 | Compact() after 50 deletions, 10⁴ faces | 4.77 ms | — | — | 3% | report | — |
| ADR-0004 | Compact() after 50 deletions, 10⁶ faces | 636 ms | — | — | 2% | report | — |
| ADR-0005 | 100 edits without Compact(), 10⁴ faces | 233 ms | — | — | 4% | report | — |
| ADR-0005 | 100 edits without Compact(), 10⁶ faces | 44478 ms | — | — | 1% | report | — |
| B0.2 | Direct call, cubic B-spline, 10⁸ evaluations | 18.4 ns | 65.4 ns | 3.6 | 4% | report | — |
| B0.2 | Direct call, circle, 10⁸ evaluations | 11.6 ns | 14.9 ns | 1.28 | 10% | report | — |
| B0.2 | Direct call, line, 10⁸ evaluations | 0.818 ns | 4.09 ns | 5.0 | 33% | report | — |
| B0.2 | VisitBuiltin against OCCT's virtual call, cubic B-spline | 22.1 ns | 62.9 ns | 2.8 | 18% | faster than OCCT | yes |
| B0.2 | VisitBuiltin against OCCT's virtual call, circle | 12.2 ns | 15 ns | 1.23 | 10% | faster than OCCT | yes |
| B0.2 | VisitBuiltin against OCCT's virtual call, line | 0.819 ns | 4.09 ns | 5.0 | 4% | faster than OCCT | yes |
| B0.2 | VisitBuiltin on a local copy, cubic B-spline | 22.1 ns | — | — | 1% | report | — |
| B0.2 | VisitBuiltin on a local copy, circle | 12.2 ns | — | — | 2% | report | — |
| B0.2 | VisitBuiltin on a local copy, line | 0.819 ns | — | — | 0% | report | — |
| B0.2 | Virtual call per evaluation, cubic B-spline | 22.5 ns | 73.6 ns | 3.3 | 15% | report | — |
| B0.2 | Virtual call per evaluation, circle | 12.2 ns | 15 ns | 1.24 | 11% | report | — |
| B0.2 | Virtual call per evaluation, line | 0.819 ns | 4.09 ns | 5.0 | 19% | report | — |
| B0.3 | OCCT's ancestor map alone, 10⁴ faces | — | 7.07 ms | — | 7% | report | — |
| B0.3 | OCCT's ancestor map alone, 10⁶ faces | — | 2707 ms | — | 3% | report | — |
| B0.3 | Edges of every vertex, 10⁴ faces | 453 µs | 13.7 ms | 30 | 8% | report | — |
| B0.3 | Edges of every vertex, 10⁶ faces | 47 ms | 3255 ms | 69 | 22% | ≥ 5× | yes |
| B0.3 | Faces of every edge, 10⁴ faces | 744 µs | 6.7 ms | 9.0 | 5% | report | — |
| B0.3 | Faces of every edge, 10⁶ faces | 103 ms | 2663 ms | 26 | 3% | ≥ 5× | yes |
| B0.3 | Build the sheet, 10⁴ faces | 6.96 ms | 163 ms | 23 | 3% | report | — |
| B0.3 | Build the sheet, 10⁶ faces | 927 ms | 16633 ms | 18 | 6% | report | — |
| B0.5 | Delete one face, 10⁴ faces | 2.31 ms | 13.2 ms | 5.7 | 14% | report | — |
| B0.5 | Delete one face, 10⁶ faces | 447 ms | 1636 ms | 3.7 | 2% | report | — |
| B0.5 | Replace one face's surface, 10⁴ faces | 2.29 ms | 13.4 ms | 5.9 | 11% | report | — |
| B0.5 | Replace one face's surface, 10⁶ faces | 439 ms | 1657 ms | 3.8 | 9% | report | — |
| B0.5 | SplitEdge/10000 | 2.3 ms | — | — | 1% | ≥ parity | no OCCT run |
| B0.5 | SplitEdge/1000000 | 438 ms | — | — | 1% | ≥ parity | no OCCT run |
| B0.5 | SplitFace/10000 | 2.31 ms | — | — | 0% | ≥ parity | no OCCT run |
| B0.5 | SplitFace/1000000 | 439 ms | — | — | 3% | ≥ parity | no OCCT run |
| B1.1 | Curve D0, cubic, 24 poles | 19.1 ns | 73.7 ns | 3.9 | 11% | ≥ parity | yes |
| B1.1 | Curve D0, rational cubic | 22.6 ns | 67.8 ns | 3.0 | 3% | ≥ parity | yes |
| B1.1 | Curve D0, degree 7 | 61.3 ns | 137 ns | 2.2 | 8% | ≥ parity | yes |
| B1.1 | Curve D0, batch of 4096, cubic, bit-identical | 21.9 ns | 7.85 ns | 0.36 | 5% | report | — |
| B1.1 | Curve D0, batch of 4096, rational cubic, bit-identical | 25.2 ns | 10.6 ns | 0.42 | 1% | report | — |
| B1.1 | Curve D0, batch of 4096, cubic, fast | 3.91 ns | 7.82 ns | 2.0 | 6% | ≥ 2× | yes |
| B1.1 | Curve D0, batch of 4096, rational cubic, fast | 5.31 ns | 10.6 ns | 1.99 | 1% | ≥ 2× | no |
| B1.1 | Curve D0–D2, cubic | 80 ns | 76.8 ns | 0.96 | 7% | ≥ parity | parity |
| B1.1 | Curve D0–D2, rational cubic | 87.2 ns | 114 ns | 1.31 | 9% | ≥ parity | yes |
| B1.1 | Surface D0, bicubic 12 × 12 | 45.8 ns | 133 ns | 2.9 | 4% | ≥ parity | yes |
| B1.1 | Surface D0, rational bicubic | 51.5 ns | 146 ns | 2.8 | 2% | ≥ parity | yes |
| B1.1 | Surface D0–D1, bicubic | 147 ns | 179 ns | 1.21 | 2% | ≥ parity | yes |
| B1.2 | Grid 256 × 256, bicubic | 9.33 ns | 33.2 ns | 3.6 | 3% | ≥ parity | yes |
| B1.2 | Grid 256 × 256, rational bicubic | 12.6 ns | 43.9 ns | 3.5 | 0% | ≥ parity | yes |
| B1.3 | Curve projection, rational cubic | 4.9 µs | 16.7 µs | 3.4 | 1% | ≥ parity | yes |
| B1.3 | Surface projection, bicubic | 92.6 µs | 99.9 µs | 1.08 | 1% | ≥ parity | parity |
| B1.4 | Approximation, 1000 points, 1e-3 mm | 1.02 ms | 23.7 ms | 23 | 1% | ≥ parity at equal tolerance | see note (not like for like) |
| B1.4 | Interpolation, 1000 points, cubic | 52 µs | 87.2 µs | 1.68 | 1% | ≥ parity | yes |
| B1.5 | Bezier decomposition, 21 segments | 3.7 µs | 5.03 µs | 1.36 | 7% | ≥ parity | yes |
| B1.5 | Degree elevation, cubic to quintic, 21 segments | 13.6 µs | 14.3 µs | 1.05 | 6% | ≥ parity | parity |
| B1.5 | 100 knot insertions | 73 µs | 175 µs | 2.4 | 1% | ≥ parity | yes |
| B1.6 | Arc length of the whole curve, rational cubic | 98.4 µs | 93 µs | 0.95 | 5% | ≥ parity at equal accuracy | parity |
| B1.6 | Frenet frame and curvature, rational cubic | 121 ns | 157 ns | 1.30 | 6% | ≥ parity | yes |
| B1.6 | Parameter at a given length | 161 µs | 244 µs | 1.51 | 12% | ≥ parity | yes |
| B1.6 | Principal curvatures and directions, bicubic | 343 ns | 308 ns | 0.90 | 4% | ≥ parity | no |

## Dispatch overhead

B0.2's second target, from the Stereon side of the rows above (ADR-0003, revised 2026-11-01): the median time per evaluation through `VisitBuiltin` over the same loop run through `VisitBuiltin` on a local copy of the geometry, which isolates what type erasure and dispatch cost. The ratio to the direct call on the concrete type is reported, not gated.

| ID | Curve | Direct | Local copy | VisitBuiltin | VisitBuiltin / copy | VisitBuiltin / direct | Target | Met |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| B0.2 | line | 0.818 ns | 0.819 ns | 0.819 ns | 1.00× | 1.00× | ≤ 1.10× | yes |
| B0.2 | circle | 11.6 ns | 12.2 ns | 12.2 ns | 1.00× | 1.06× | ≤ 1.10× | yes |
| B0.2 | cubic B-spline | 18.4 ns | 22.1 ns | 22.1 ns | 1.00× | 1.20× | ≤ 1.10× | yes |

## Memory and copying

Median over the repetitions of each counter. Live heap is the bytes the sheet holds in live allocations with the shape alive, measured the same way for both libraries (bench/Memory.hpp); pages copied are the pages of the edited snapshot not shared with the body (ADR-0004, ADR-0005). Stereon / OCCT below 100 % means Stereon uses less.

| ID | Measure | Stereon | OCCT | Stereon / OCCT | Target | Met |
| --- | --- | --- | --- | --- | --- | --- |
| B0.4 | Live heap per face, 10⁴ faces | 617 B | 1,611 B | 38% | report | — |
| B0.4 | Live heap per edge, 10⁴ faces | 305 B | 797 B | 38% | report | — |
| B0.4 | Table pages per face, 10⁴ faces | 613 B | — | — | report | — |
| B0.5 | Pages copied, one surface replaced, 10⁴ faces | 0.361% | — | — | report | — |
| B0.5 | Heap added, one surface replaced, 10⁴ faces | 38,144 B | 479,952 B | 8% | report | — |
| B0.5 | Pages copied, one face deleted, 10⁴ faces | 0.963% | — | — | report | — |
| B0.5 | Heap added, one face deleted, 10⁴ faces | 102,816 B | 479,856 B | 21% | report | — |
| ADR-0005 | Body and its edited copy over the body, 100 edits, 10⁴ faces | 1.333× | — | — | report | — |
| ADR-0004 | Compact() result over a fresh body, after 50 deletions, 10⁴ faces | 0.995× | — | — | report | — |
| B0.4 | Live heap per face, 10⁶ faces | 608 B | 1,601 B | 38% | < 50% of OCCT | yes |
| B0.4 | Live heap per edge, 10⁶ faces | 304 B | 800 B | 38% | report | — |
| B0.4 | Table pages per face, 10⁶ faces | 604 B | — | — | report | — |
| B0.4 |   of which vertices | 40 B | — | — | report | — |
| B0.4 |   of which edges | 144 B | — | — | report | — |
| B0.4 |   of which coedges | 128 B | — | — | report | — |
| B0.4 |   of which loops and faces | 44 B | — | — | report | — |
| B0.4 |   of which curve pool | 104 B | — | — | report | — |
| B0.4 |   of which pcurve pool | 144 B | — | — | report | — |
| B0.4 |   of which surface pool | 0 B | — | — | report | — |
| B0.5 | Pages copied, one surface replaced, 10⁶ faces | 0.004% | — | — | < 1% of body | yes |
| B0.5 | Heap added, one surface replaced, 10⁶ faces | 1,339,936 B | 48,000,432 B | 3% | report | — |
| B0.5 | Pages copied, one face deleted, 10⁶ faces | 0.010% | — | — | < 1% of body | yes |
| B0.5 | Heap added, one face deleted, 10⁶ faces | 1,404,640 B | 48,000,048 B | 3% | report | — |
| ADR-0005 | Body and its edited copy over the body, 100 edits, 10⁶ faces | 1.005× | — | — | < 1.2× | yes |
| ADR-0004 | Compact() result over a fresh body, after 50 deletions, 10⁶ faces | 1.000× | — | — | < 1.05× | yes |

## Notes

- **B0.2/Direct/bspline**, `sum_x_mm`: 229372.50000000178 (Stereon), 229372.50000000175 (OCCT).
- **B0.2/Direct/circle**, `sum_x_mm`: 65535.999999999927 (OCCT), 65535.999999999927 (Stereon).
- **B0.2/Direct/line**, `sum_x_mm`: 168478.13727650364 (OCCT), 168478.13727650364 (Stereon).
- **B0.2/Dispatch/bspline**, `sum_x_mm`: 229372.50000000178 (Stereon), 229372.50000000175 (OCCT).
- **B0.2/Dispatch/circle**, `sum_x_mm`: 65535.999999999927 (OCCT), 65535.999999999927 (Stereon).
- **B0.2/Dispatch/line**, `sum_x_mm`: 168478.13727650364 (OCCT), 168478.13727650364 (Stereon).
- **B0.2/DispatchCopy/bspline**, `sum_x_mm`: 229372.50000000178 (Stereon).
- **B0.2/DispatchCopy/circle**, `sum_x_mm`: 65535.999999999927 (Stereon).
- **B0.2/DispatchCopy/line**, `sum_x_mm`: 168478.13727650364 (Stereon).
- **B0.2/Virtual/bspline**, `sum_x_mm`: 229372.50000000178 (Stereon), 229372.50000000175 (OCCT).
- **B0.2/Virtual/circle**, `sum_x_mm`: 65535.999999999927 (OCCT), 65535.999999999927 (Stereon).
- **B0.2/Virtual/line**, `sum_x_mm`: 168478.13727650364 (Stereon), 168478.13727650364 (OCCT).
- **B1.6/ArcLength/rational-cubic**, `length_mm`: 111.80161013322777 (Stereon), 111.80161013322775 (OCCT).
- **Approximation is not like for like.** Stereon returns the fewest-pole least-squares cubic within the budget; OCCT's `GeomAPI_PointsToBSpline` also optimises parameters and continuity. Read the ratio as indicative.
