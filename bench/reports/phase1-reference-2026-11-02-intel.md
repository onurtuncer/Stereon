# Benchmark reference run — dispatch (B0.2), topology (B0.3–B0.5), curves and surfaces (B1.1–B1.6)

## Setup

| | |
| --- | --- |
| Date | 2026-11-02 |
| Commit | ab41d1a124182bbe874daafdaeaecef6bef3e68f |
| Ref | main |
| Workflow run | https://github.com/onurtuncer/Stereon/actions/runs/41434476419 |
| Runner | GitHub Actions ubuntu-24.04 (Linux, X64), container gcc:16.2 |
| CPU | INTEL(R) XEON(R) PLATINUM 8573C |
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
| ADR-0004 | Compact() after 50 deletions, 10⁴ faces | 5.57 ms | — | — | 0% | report | — |
| ADR-0004 | Compact() after 50 deletions, 10⁶ faces | 843 ms | — | — | 1% | report | — |
| ADR-0005 | 100 edits without Compact(), 10⁴ faces | 285 ms | — | — | 1% | report | — |
| ADR-0005 | 100 edits without Compact(), 10⁶ faces | 49689 ms | — | — | 0% | report | — |
| B0.2 | Direct call, cubic B-spline, 10⁸ evaluations | 19 ns | 68.4 ns | 3.6 | 0% | report | — |
| B0.2 | Direct call, circle, 10⁸ evaluations | 10.6 ns | 13 ns | 1.23 | 1% | report | — |
| B0.2 | Direct call, line, 10⁸ evaluations | 0.673 ns | 4.04 ns | 6.0 | 0% | report | — |
| B0.2 | VisitBuiltin against OCCT's virtual call, cubic B-spline | 22.8 ns | 69 ns | 3.0 | 1% | faster than OCCT | yes |
| B0.2 | VisitBuiltin against OCCT's virtual call, circle | 10.9 ns | 13 ns | 1.19 | 1% | faster than OCCT | yes |
| B0.2 | VisitBuiltin against OCCT's virtual call, line | 0.673 ns | 4.04 ns | 6.0 | 1% | faster than OCCT | yes |
| B0.2 | VisitBuiltin on a local copy, cubic B-spline | 22.8 ns | — | — | 1% | report | — |
| B0.2 | VisitBuiltin on a local copy, circle | 10.9 ns | — | — | 1% | report | — |
| B0.2 | VisitBuiltin on a local copy, line | 0.673 ns | — | — | 2% | report | — |
| B0.2 | Virtual call per evaluation, cubic B-spline | 23.2 ns | 69 ns | 3.0 | 1% | report | — |
| B0.2 | Virtual call per evaluation, circle | 11.2 ns | 13 ns | 1.16 | 1% | report | — |
| B0.2 | Virtual call per evaluation, line | 0.822 ns | 4.04 ns | 4.9 | 1% | report | — |
| B0.3 | OCCT's ancestor map alone, 10⁴ faces | — | 8.94 ms | — | 5% | report | — |
| B0.3 | OCCT's ancestor map alone, 10⁶ faces | — | 2727 ms | — | 3% | report | — |
| B0.3 | Edges of every vertex, 10⁴ faces | 694 µs | 16.4 ms | 24 | 5% | report | — |
| B0.3 | Edges of every vertex, 10⁶ faces | 79.2 ms | 4175 ms | 53 | 1% | ≥ 5× | yes |
| B0.3 | Faces of every edge, 10⁴ faces | 896 µs | 9.06 ms | 10 | 3% | report | — |
| B0.3 | Faces of every edge, 10⁶ faces | 114 ms | 2739 ms | 24 | 26% | ≥ 5× | yes |
| B0.3 | Build the sheet, 10⁴ faces | 8.49 ms | 262 ms | 31 | 3% | report | — |
| B0.3 | Build the sheet, 10⁶ faces | 1119 ms | 25536 ms | 23 | 23% | report | — |
| B0.5 | Delete one face, 10⁴ faces | 2.83 ms | 29.8 ms | 11 | 1% | report | — |
| B0.5 | Delete one face, 10⁶ faces | 498 ms | 3072 ms | 6.2 | 1% | report | — |
| B0.5 | Replace one face's surface, 10⁴ faces | 2.79 ms | 29.8 ms | 11 | 2% | report | — |
| B0.5 | Replace one face's surface, 10⁶ faces | 490 ms | 3089 ms | 6.3 | 1% | report | — |
| B0.5 | SplitEdge/10000 | 2.8 ms | — | — | 1% | ≥ parity | no OCCT run |
| B0.5 | SplitEdge/1000000 | 489 ms | — | — | 1% | ≥ parity | no OCCT run |
| B0.5 | SplitFace/10000 | 2.83 ms | — | — | 2% | ≥ parity | no OCCT run |
| B0.5 | SplitFace/1000000 | 488 ms | — | — | 1% | ≥ parity | no OCCT run |
| B1.1 | Curve D0, cubic, 24 poles | 19.6 ns | 70.6 ns | 3.6 | 2% | ≥ parity | yes |
| B1.1 | Curve D0, rational cubic | 26.1 ns | 76.2 ns | 2.9 | 1% | ≥ parity | yes |
| B1.1 | Curve D0, degree 7 | 59.4 ns | 132 ns | 2.2 | 1% | ≥ parity | yes |
| B1.1 | Curve D0, batch of 4096, cubic, bit-identical | 17 ns | 8.76 ns | 0.52 | 2% | report | — |
| B1.1 | Curve D0, batch of 4096, rational cubic, bit-identical | 18.8 ns | 10.4 ns | 0.55 | 1% | report | — |
| B1.1 | Curve D0, batch of 4096, cubic, fast | 4.88 ns | 8.76 ns | 1.79 | 1% | ≥ 2× | no |
| B1.1 | Curve D0, batch of 4096, rational cubic, fast | 5.89 ns | 10.4 ns | 1.77 | 1% | ≥ 2× | no |
| B1.1 | Curve D0–D2, cubic | 70 ns | 79.3 ns | 1.13 | 1% | ≥ parity | yes |
| B1.1 | Curve D0–D2, rational cubic | 78.2 ns | 117 ns | 1.50 | 1% | ≥ parity | yes |
| B1.1 | Surface D0, bicubic 12 × 12 | 54.2 ns | 135 ns | 2.5 | 1% | ≥ parity | yes |
| B1.1 | Surface D0, rational bicubic | 61.7 ns | 166 ns | 2.7 | 1% | ≥ parity | yes |
| B1.1 | Surface D0–D1, bicubic | 146 ns | 213 ns | 1.47 | 1% | ≥ parity | yes |
| B1.2 | Grid 256 × 256, bicubic | 12.5 ns | 39.8 ns | 3.2 | 1% | ≥ parity | yes |
| B1.2 | Grid 256 × 256, rational bicubic | 15.8 ns | 45.7 ns | 2.9 | 2% | ≥ parity | yes |
| B1.3 | Curve projection, rational cubic | 6.02 µs | 18.5 µs | 3.1 | 8% | ≥ parity | yes |
| B1.3 | Surface projection, bicubic | 105 µs | 116 µs | 1.10 | 1% | ≥ parity | yes |
| B1.4 | Approximation, 1000 points, 1e-3 mm | 1.2 ms | 32.1 ms | 27 | 1% | ≥ parity at equal tolerance | see note (not like for like) |
| B1.4 | Interpolation, 1000 points, cubic | 75.3 µs | 115 µs | 1.52 | 2% | ≥ parity | yes |
| B1.5 | Bezier decomposition, 21 segments | 4.84 µs | 6.61 µs | 1.37 | 4% | ≥ parity | yes |
| B1.5 | Degree elevation, cubic to quintic, 21 segments | 21.3 µs | 16.5 µs | 0.78 | 2% | ≥ parity | no |
| B1.5 | 100 knot insertions | 86.8 µs | 230 µs | 2.7 | 3% | ≥ parity | yes |
| B1.6 | Arc length of the whole curve, rational cubic | 83.1 µs | 99.3 µs | 1.20 | 1% | ≥ parity at equal accuracy | yes |
| B1.6 | Frenet frame and curvature, rational cubic | 115 ns | 162 ns | 1.41 | 7% | ≥ parity | yes |
| B1.6 | Parameter at a given length | 137 µs | 258 µs | 1.89 | 1% | ≥ parity | yes |
| B1.6 | Principal curvatures and directions, bicubic | 371 ns | 366 ns | 0.99 | 1% | ≥ parity | parity |

## Dispatch overhead

B0.2's second target, from the Stereon side of the rows above (ADR-0003, revised 2026-11-01): the median time per evaluation through `VisitBuiltin` over the same loop run through `VisitBuiltin` on a local copy of the geometry, which isolates what type erasure and dispatch cost. The ratio to the direct call on the concrete type is reported, not gated.

| ID | Curve | Direct | Local copy | VisitBuiltin | VisitBuiltin / copy | VisitBuiltin / direct | Target | Met |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| B0.2 | line | 0.673 ns | 0.673 ns | 0.673 ns | 1.00× | 1.00× | ≤ 1.10× | yes |
| B0.2 | circle | 10.6 ns | 10.9 ns | 10.9 ns | 1.00× | 1.03× | ≤ 1.10× | yes |
| B0.2 | cubic B-spline | 19 ns | 22.8 ns | 22.8 ns | 1.00× | 1.20× | ≤ 1.10× | yes |

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
- **B0.2/Direct/circle**, `sum_x_mm`: 65535.999999999927 (Stereon), 65535.999999999927 (OCCT).
- **B0.2/Direct/line**, `sum_x_mm`: 168478.13727650364 (Stereon), 168478.13727650364 (OCCT).
- **B0.2/Dispatch/bspline**, `sum_x_mm`: 229372.50000000175 (OCCT), 229372.50000000178 (Stereon).
- **B0.2/Dispatch/circle**, `sum_x_mm`: 65535.999999999927 (OCCT), 65535.999999999927 (Stereon).
- **B0.2/Dispatch/line**, `sum_x_mm`: 168478.13727650364 (Stereon), 168478.13727650364 (OCCT).
- **B0.2/DispatchCopy/bspline**, `sum_x_mm`: 229372.50000000178 (Stereon).
- **B0.2/DispatchCopy/circle**, `sum_x_mm`: 65535.999999999927 (Stereon).
- **B0.2/DispatchCopy/line**, `sum_x_mm`: 168478.13727650364 (Stereon).
- **B0.2/Virtual/bspline**, `sum_x_mm`: 229372.50000000178 (Stereon), 229372.50000000175 (OCCT).
- **B0.2/Virtual/circle**, `sum_x_mm`: 65535.999999999927 (OCCT), 65535.999999999927 (Stereon).
- **B0.2/Virtual/line**, `sum_x_mm`: 168478.13727650364 (OCCT), 168478.13727650364 (Stereon).
- **B1.6/ArcLength/rational-cubic**, `length_mm`: 111.80161013322775 (OCCT), 111.80161013322777 (Stereon).
- **Approximation is not like for like.** Stereon returns the fewest-pole least-squares cubic within the budget; OCCT's `GeomAPI_PointsToBSpline` also optimises parameters and continuity. Read the ratio as indicative.
