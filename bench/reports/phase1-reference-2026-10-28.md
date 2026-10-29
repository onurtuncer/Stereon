# Benchmark reference run — dispatch (B0.2), topology (B0.3–B0.5), curves and surfaces (B1.1–B1.6)

## Setup

| | |
| --- | --- |
| Date | 2026-10-28 |
| Commit | 8e3aa2ef224724d385f5d91a37defce89c277622 |
| Ref | 149/merge |
| Workflow run | https://github.com/onurtuncer/Stereon/actions/runs/40738691891 |
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
| ADR-0004 | Compact() after 50 deletions, 10⁴ faces | 5.85 ms | — | — | 49% | report | — |
| ADR-0004 | Compact() after 50 deletions, 10⁶ faces | 719 ms | — | — | 40% | report | — |
| ADR-0005 | 100 edits without Compact(), 10⁴ faces | 293 ms | — | — | 3% | report | — |
| ADR-0005 | 100 edits without Compact(), 10⁶ faces | 44369 ms | — | — | 3% | report | — |
| B0.2 | Direct call, cubic B-spline, 10⁸ evaluations | 50.9 ns | 81.2 ns | 1.59 | 1% | report | — |
| B0.2 | Direct call, circle, 10⁸ evaluations | 14.8 ns | 18.7 ns | 1.26 | 1% | report | — |
| B0.2 | Direct call, line, 10⁸ evaluations | 0.938 ns | 5 ns | 5.3 | 12% | report | — |
| B0.2 | VisitBuiltin against OCCT's virtual call, cubic B-spline | 55.7 ns | 80.8 ns | 1.45 | 0% | faster than OCCT | yes |
| B0.2 | VisitBuiltin against OCCT's virtual call, circle | 15 ns | 18.9 ns | 1.26 | 2% | faster than OCCT | yes |
| B0.2 | VisitBuiltin against OCCT's virtual call, line | 0.938 ns | 4.69 ns | 5.0 | 0% | faster than OCCT | yes |
| B0.2 | Virtual call per evaluation, cubic B-spline | 55.8 ns | 80.8 ns | 1.45 | 0% | report | — |
| B0.2 | Virtual call per evaluation, circle | 15.3 ns | 18.9 ns | 1.24 | 3% | report | — |
| B0.2 | Virtual call per evaluation, line | 0.941 ns | 4.69 ns | 5.0 | 2% | report | — |
| B0.3 | OCCT's ancestor map alone, 10⁴ faces | — | 7.99 ms | — | 8% | report | — |
| B0.3 | OCCT's ancestor map alone, 10⁶ faces | — | 2643 ms | — | 2% | report | — |
| B0.3 | Edges of every vertex, 10⁴ faces | 530 µs | 13.9 ms | 26 | 57% | report | — |
| B0.3 | Edges of every vertex, 10⁶ faces | 53.9 ms | 3288 ms | 61 | 21% | ≥ 5× | yes |
| B0.3 | Faces of every edge, 10⁴ faces | 953 µs | 7.87 ms | 8.3 | 3% | report | — |
| B0.3 | Faces of every edge, 10⁶ faces | 114 ms | 2617 ms | 23 | 1% | ≥ 5× | yes |
| B0.3 | Build the sheet, 10⁴ faces | 8.9 ms | 223 ms | 25 | 2% | report | — |
| B0.3 | Build the sheet, 10⁶ faces | 1088 ms | 22297 ms | 20 | 4% | report | — |
| B0.5 | Delete one face, 10⁴ faces | 2.86 ms | 16.3 ms | 5.7 | 17% | report | — |
| B0.5 | Delete one face, 10⁶ faces | 440 ms | 1953 ms | 4.4 | 2% | report | — |
| B0.5 | Replace one face's surface, 10⁴ faces | 2.84 ms | 17 ms | 6.0 | 15% | report | — |
| B0.5 | Replace one face's surface, 10⁶ faces | 435 ms | 1970 ms | 4.5 | 4% | report | — |
| B0.5 | SplitEdge/10000 | 2.85 ms | — | — | 1% | ≥ parity | no OCCT run |
| B0.5 | SplitEdge/1000000 | 435 ms | — | — | 1% | ≥ parity | no OCCT run |
| B0.5 | SplitFace/10000 | 2.87 ms | — | — | 0% | ≥ parity | no OCCT run |
| B0.5 | SplitFace/1000000 | 437 ms | — | — | 1% | ≥ parity | no OCCT run |
| B1.1 | Curve D0, cubic, 24 poles | 51.5 ns | 84.9 ns | 1.65 | 1% | ≥ parity | yes |
| B1.1 | Curve D0, rational cubic | 57.3 ns | 88.4 ns | 1.54 | 1% | ≥ parity | yes |
| B1.1 | Curve D0, degree 7 | 94.7 ns | 153 ns | 1.61 | 1% | ≥ parity | yes |
| B1.1 | Curve D0, batch of 4096, cubic, bit-identical | 25.9 ns | 11 ns | 0.43 | 3% | report | — |
| B1.1 | Curve D0, batch of 4096, rational cubic, bit-identical | 29.6 ns | 15 ns | 0.51 | 1% | report | — |
| B1.1 | Curve D0, batch of 4096, cubic, fast | 4.23 ns | 11 ns | 2.6 | 1% | ≥ 2× | yes |
| B1.1 | Curve D0, batch of 4096, rational cubic, fast | 6.02 ns | 15 ns | 2.5 | 3% | ≥ 2× | yes |
| B1.1 | Curve D0–D2, cubic | 94.4 ns | 99.7 ns | 1.06 | 1% | ≥ parity | parity |
| B1.1 | Curve D0–D2, rational cubic | 105 ns | 139 ns | 1.33 | 2% | ≥ parity | yes |
| B1.1 | Surface D0, bicubic 12 × 12 | 109 ns | 173 ns | 1.58 | 14% | ≥ parity | yes |
| B1.1 | Surface D0, rational bicubic | 115 ns | 198 ns | 1.73 | 0% | ≥ parity | yes |
| B1.1 | Surface D0–D1, bicubic | 191 ns | 255 ns | 1.33 | 3% | ≥ parity | yes |
| B1.2 | Grid 256 × 256, bicubic | 12.4 ns | 42 ns | 3.4 | 1% | ≥ parity | yes |
| B1.2 | Grid 256 × 256, rational bicubic | 17.1 ns | 56.4 ns | 3.3 | 1% | ≥ parity | yes |
| B1.3 | Curve projection, rational cubic | 5.51 µs | 20.6 µs | 3.7 | 1% | ≥ parity | yes |
| B1.3 | Surface projection, bicubic | 104 µs | 129 µs | 1.24 | 1% | ≥ parity | yes |
| B1.4 | Approximation, 1000 points, 1e-3 mm | 1.95 ms | 30 ms | 15 | 2% | ≥ parity at equal tolerance | see note (not like for like) |
| B1.4 | Interpolation, 1000 points, cubic | 84.4 µs | 113 µs | 1.34 | 2% | ≥ parity | yes |
| B1.5 | Bezier decomposition, 21 segments | 4.99 µs | 7.21 µs | 1.44 | 11% | ≥ parity | yes |
| B1.5 | Degree elevation, cubic to quintic, 21 segments | 5.96 µs | 20.2 µs | 3.4 | 6% | ≥ parity | yes |
| B1.5 | 100 knot insertions | 90.7 µs | 229 µs | 2.5 | 2% | ≥ parity | yes |
| B1.6 | Arc length of the whole curve, rational cubic | 107 µs | 113 µs | 1.05 | 2% | ≥ parity at equal accuracy | parity |
| B1.6 | Frenet frame and curvature, rational cubic | 132 ns | 188 ns | 1.42 | 0% | ≥ parity | yes |
| B1.6 | Parameter at a given length | 186 µs | 294 µs | 1.58 | 1% | ≥ parity | yes |
| B1.6 | Principal curvatures and directions, bicubic | 444 ns | 415 ns | 0.93 | 7% | ≥ parity | parity |

## Dispatch overhead

B0.2's second target, from the Stereon side of two rows above: the median time per evaluation through `VisitBuiltin` over the direct call on the concrete type (ADR-0003).

| ID | Curve | Direct | VisitBuiltin | VisitBuiltin / direct | Target | Met |
| --- | --- | --- | --- | --- | --- | --- |
| B0.2 | line | 0.938 ns | 0.938 ns | 1.00× | ≤ 1.10× | yes |
| B0.2 | circle | 14.8 ns | 15 ns | 1.01× | ≤ 1.10× | yes |
| B0.2 | cubic B-spline | 50.9 ns | 55.7 ns | 1.09× | ≤ 1.10× | yes |

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

- **B0.2/Direct/bspline**, `sum_x_mm`: 229372.50000000175 (OCCT), 229372.50000000178 (Stereon).
- **B0.2/Direct/circle**, `sum_x_mm`: 65535.999999999927 (OCCT), 65535.999999999927 (Stereon).
- **B0.2/Direct/line**, `sum_x_mm`: 168478.13727650364 (OCCT), 168478.13727650364 (Stereon).
- **B0.2/Dispatch/bspline**, `sum_x_mm`: 229372.50000000178 (Stereon), 229372.50000000175 (OCCT).
- **B0.2/Dispatch/circle**, `sum_x_mm`: 65535.999999999927 (Stereon), 65535.999999999927 (OCCT).
- **B0.2/Dispatch/line**, `sum_x_mm`: 168478.13727650364 (Stereon), 168478.13727650364 (OCCT).
- **B0.2/Virtual/bspline**, `sum_x_mm`: 229372.50000000175 (OCCT), 229372.50000000178 (Stereon).
- **B0.2/Virtual/circle**, `sum_x_mm`: 65535.999999999927 (OCCT), 65535.999999999927 (Stereon).
- **B0.2/Virtual/line**, `sum_x_mm`: 168478.13727650364 (Stereon), 168478.13727650364 (OCCT).
- **B1.6/ArcLength/rational-cubic**, `length_mm`: 111.80161013322777 (Stereon), 111.80161013322775 (OCCT).
- **Approximation is not like for like.** Stereon returns the fewest-pole least-squares cubic within the budget; OCCT's `GeomAPI_PointsToBSpline` also optimises parameters and continuity. Read the ratio as indicative.
