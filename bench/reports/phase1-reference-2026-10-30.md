# Benchmark reference run — dispatch (B0.2), topology (B0.3–B0.5), curves and surfaces (B1.1–B1.6)

## Setup

| | |
| --- | --- |
| Date | 2026-10-30 |
| Commit | fb7366ae74dc8864d09d0f21e2632f6d1cd573b5 |
| Ref | main |
| Workflow run | https://github.com/onurtuncer/Stereon/actions/runs/41115492314 |
| Runner | GitHub Actions ubuntu-24.04 (Linux, X64), container gcc:16.2 |
| CPU | AMD EPYC 7763 64-Core Processor |
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
| ADR-0004 | Compact() after 50 deletions, 10⁴ faces | 5.81 ms | — | — | 2% | report | — |
| ADR-0004 | Compact() after 50 deletions, 10⁶ faces | 731 ms | — | — | 2% | report | — |
| ADR-0005 | 100 edits without Compact(), 10⁴ faces | 293 ms | — | — | 1% | report | — |
| ADR-0005 | 100 edits without Compact(), 10⁶ faces | 44862 ms | — | — | 3% | report | — |
| B0.2 | Direct call, cubic B-spline, 10⁸ evaluations | 22.8 ns | 80.7 ns | 3.5 | 1% | report | — |
| B0.2 | Direct call, circle, 10⁸ evaluations | 14.9 ns | 18.7 ns | 1.25 | 1% | report | — |
| B0.2 | Direct call, line, 10⁸ evaluations | 0.938 ns | 5 ns | 5.3 | 6% | report | — |
| B0.2 | VisitBuiltin against OCCT's virtual call, cubic B-spline | 27.7 ns | 81.1 ns | 2.9 | 3% | faster than OCCT | yes |
| B0.2 | VisitBuiltin against OCCT's virtual call, circle | 15 ns | 18.9 ns | 1.26 | 1% | faster than OCCT | yes |
| B0.2 | VisitBuiltin against OCCT's virtual call, line | 0.939 ns | 4.69 ns | 5.0 | 0% | faster than OCCT | yes |
| B0.2 | Virtual call per evaluation, cubic B-spline | 28.5 ns | 81.1 ns | 2.8 | 2% | report | — |
| B0.2 | Virtual call per evaluation, circle | 15.3 ns | 18.9 ns | 1.23 | 1% | report | — |
| B0.2 | Virtual call per evaluation, line | 0.941 ns | 4.68 ns | 5.0 | 7% | report | — |
| B0.3 | OCCT's ancestor map alone, 10⁴ faces | — | 7.71 ms | — | 9% | report | — |
| B0.3 | OCCT's ancestor map alone, 10⁶ faces | — | 2565 ms | — | 22% | report | — |
| B0.3 | Edges of every vertex, 10⁴ faces | 526 µs | 13.9 ms | 26 | 24% | report | — |
| B0.3 | Edges of every vertex, 10⁶ faces | 53.3 ms | 3385 ms | 64 | 3% | ≥ 5× | yes |
| B0.3 | Faces of every edge, 10⁴ faces | 1.04 ms | 7.8 ms | 7.5 | 18% | report | — |
| B0.3 | Faces of every edge, 10⁶ faces | 125 ms | 2603 ms | 21 | 5% | ≥ 5× | yes |
| B0.3 | Build the sheet, 10⁴ faces | 8.86 ms | 227 ms | 26 | 2% | report | — |
| B0.3 | Build the sheet, 10⁶ faces | 1073 ms | 22931 ms | 21 | 1% | report | — |
| B0.5 | Delete one face, 10⁴ faces | 2.87 ms | 16.1 ms | 5.6 | 2% | report | — |
| B0.5 | Delete one face, 10⁶ faces | 444 ms | 1955 ms | 4.4 | 4% | report | — |
| B0.5 | Replace one face's surface, 10⁴ faces | 2.83 ms | 16.1 ms | 5.7 | 13% | report | — |
| B0.5 | Replace one face's surface, 10⁶ faces | 444 ms | 1953 ms | 4.4 | 2% | report | — |
| B0.5 | SplitEdge/10000 | 2.83 ms | — | — | 1% | ≥ parity | no OCCT run |
| B0.5 | SplitEdge/1000000 | 441 ms | — | — | 3% | ≥ parity | no OCCT run |
| B0.5 | SplitFace/10000 | 2.87 ms | — | — | 1% | ≥ parity | no OCCT run |
| B0.5 | SplitFace/1000000 | 445 ms | — | — | 2% | ≥ parity | no OCCT run |
| B1.1 | Curve D0, cubic, 24 poles | 23.8 ns | 84.8 ns | 3.6 | 0% | ≥ parity | yes |
| B1.1 | Curve D0, rational cubic | 27.6 ns | 88.3 ns | 3.2 | 0% | ≥ parity | yes |
| B1.1 | Curve D0, degree 7 | 78.5 ns | 152 ns | 1.94 | 1% | ≥ parity | yes |
| B1.1 | Curve D0, batch of 4096, cubic, bit-identical | 26.5 ns | 11 ns | 0.42 | 3% | report | — |
| B1.1 | Curve D0, batch of 4096, rational cubic, bit-identical | 30.1 ns | 15 ns | 0.50 | 3% | report | — |
| B1.1 | Curve D0, batch of 4096, cubic, fast | 4.82 ns | 11 ns | 2.3 | 1% | ≥ 2× | yes |
| B1.1 | Curve D0, batch of 4096, rational cubic, fast | 6.43 ns | 15 ns | 2.3 | 0% | ≥ 2× | yes |
| B1.1 | Curve D0–D2, cubic | 97.1 ns | 98.6 ns | 1.02 | 1% | ≥ parity | parity |
| B1.1 | Curve D0–D2, rational cubic | 105 ns | 139 ns | 1.32 | 1% | ≥ parity | yes |
| B1.1 | Surface D0, bicubic 12 × 12 | 62.5 ns | 172 ns | 2.8 | 0% | ≥ parity | yes |
| B1.1 | Surface D0, rational bicubic | 68.3 ns | 198 ns | 2.9 | 0% | ≥ parity | yes |
| B1.1 | Surface D0–D1, bicubic | 191 ns | 255 ns | 1.33 | 3% | ≥ parity | yes |
| B1.2 | Grid 256 × 256, bicubic | 12.6 ns | 42.6 ns | 3.4 | 2% | ≥ parity | yes |
| B1.2 | Grid 256 × 256, rational bicubic | 16.5 ns | 56.4 ns | 3.4 | 1% | ≥ parity | yes |
| B1.3 | Curve projection, rational cubic | 5.67 µs | 20.6 µs | 3.6 | 0% | ≥ parity | yes |
| B1.3 | Surface projection, bicubic | 104 µs | 129 µs | 1.24 | 1% | ≥ parity | yes |
| B1.4 | Approximation, 1000 points, 1e-3 mm | 1.27 ms | 29.9 ms | 24 | 12% | ≥ parity at equal tolerance | see note (not like for like) |
| B1.4 | Interpolation, 1000 points, cubic | 63.8 µs | 113 µs | 1.77 | 2% | ≥ parity | yes |
| B1.5 | Bezier decomposition, 21 segments | 4.96 µs | 7.19 µs | 1.45 | 9% | ≥ parity | yes |
| B1.5 | Degree elevation, cubic to quintic, 21 segments | 15.8 µs | 20.2 µs | 1.28 | 3% | ≥ parity | yes |
| B1.5 | 100 knot insertions | 96.3 µs | 229 µs | 2.4 | 3% | ≥ parity | yes |
| B1.6 | Arc length of the whole curve, rational cubic | 108 µs | 114 µs | 1.06 | 1% | ≥ parity at equal accuracy | parity |
| B1.6 | Frenet frame and curvature, rational cubic | 134 ns | 187 ns | 1.40 | 0% | ≥ parity | yes |
| B1.6 | Parameter at a given length | 187 µs | 293 µs | 1.57 | 2% | ≥ parity | yes |
| B1.6 | Principal curvatures and directions, bicubic | 443 ns | 398 ns | 0.90 | 0% | ≥ parity | no |

## Dispatch overhead

B0.2's second target, from the Stereon side of two rows above: the median time per evaluation through `VisitBuiltin` over the direct call on the concrete type (ADR-0003).

| ID | Curve | Direct | VisitBuiltin | VisitBuiltin / direct | Target | Met |
| --- | --- | --- | --- | --- | --- | --- |
| B0.2 | line | 0.938 ns | 0.939 ns | 1.00× | ≤ 1.10× | yes |
| B0.2 | circle | 14.9 ns | 15 ns | 1.01× | ≤ 1.10× | yes |
| B0.2 | cubic B-spline | 22.8 ns | 27.7 ns | 1.22× | ≤ 1.10× | no |

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
| B0.5 | Heap added, one face deleted, 10⁶ faces | 1,404,608 B | 48,000,048 B | 3% | report | — |
| ADR-0005 | Body and its edited copy over the body, 100 edits, 10⁶ faces | 1.005× | — | — | < 1.2× | yes |
| ADR-0004 | Compact() result over a fresh body, after 50 deletions, 10⁶ faces | 1.000× | — | — | < 1.05× | yes |

## Notes

- **B0.2/Direct/bspline**, `sum_x_mm`: 229372.50000000178 (Stereon), 229372.50000000175 (OCCT).
- **B0.2/Direct/circle**, `sum_x_mm`: 65535.999999999927 (Stereon), 65535.999999999927 (OCCT).
- **B0.2/Direct/line**, `sum_x_mm`: 168478.13727650364 (OCCT), 168478.13727650364 (Stereon).
- **B0.2/Dispatch/bspline**, `sum_x_mm`: 229372.50000000178 (Stereon), 229372.50000000175 (OCCT).
- **B0.2/Dispatch/circle**, `sum_x_mm`: 65535.999999999927 (Stereon), 65535.999999999927 (OCCT).
- **B0.2/Dispatch/line**, `sum_x_mm`: 168478.13727650364 (OCCT), 168478.13727650364 (Stereon).
- **B0.2/Virtual/bspline**, `sum_x_mm`: 229372.50000000178 (Stereon), 229372.50000000175 (OCCT).
- **B0.2/Virtual/circle**, `sum_x_mm`: 65535.999999999927 (Stereon), 65535.999999999927 (OCCT).
- **B0.2/Virtual/line**, `sum_x_mm`: 168478.13727650364 (OCCT), 168478.13727650364 (Stereon).
- **B1.6/ArcLength/rational-cubic**, `length_mm`: 111.80161013322775 (OCCT), 111.80161013322777 (Stereon).
- **Approximation is not like for like.** Stereon returns the fewest-pole least-squares cubic within the budget; OCCT's `GeomAPI_PointsToBSpline` also optimises parameters and continuity. Read the ratio as indicative.
