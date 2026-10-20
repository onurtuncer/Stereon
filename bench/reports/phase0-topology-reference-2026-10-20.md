# Benchmark reference run — topology (B0.3–B0.5), curves and surfaces (B1.1–B1.6)

## Setup

| | |
| --- | --- |
| Date | 2026-10-20 |
| Commit | 2c24d79099d14f83b44e1eb99034f9e0fca7bdd4 |
| Ref | 77/merge |
| Workflow run | https://github.com/onurtuncer/Stereon/actions/runs/39399395722 |
| Runner | GitHub Actions ubuntu-24.04 (Linux, X64), container gcc:16.2 |
| CPU | AMD EPYC 7763 64-Core Processor |
| CPUs | 4 |
| Memory | 15.6 GB |
| Kernel | Linux 6.17.0-1022-azure |
| Compiler | g++ (GCC) 16.2.0 |
| OCCT tag | V8_0_1 |
| OCCT commit | b8f597c677811d1f9f4d8a97f5ae2825c0353a42 |
| OCCT compiler | c++ (GCC) 16.2.0 |
| OCCT build type | Release |
| Repetitions | 5, 15 |
| memory_measure | live heap, glibc mallinfo2 (bench/Memory.hpp) |

## Results

Median CPU time over the repetitions, both libraries interleaved: per evaluated parameter, projected point, property or length query; per fit (B1.4), per 100 insertions and per decomposition (B1.5); per build, per query over the whole sheet or per edit (B0.3, B0.5, ADR-0004, ADR-0005). Ratio = OCCT / Stereon, so > 1 means Stereon is faster; within 10 % of 1 is parity, and a row with a ratio target of its own (B0.3: ≥ 5×) is met when it reaches it. A target of "report" is not a gate. Spread is the largest relative deviation of any repetition from its median, over both libraries.

| ID | Benchmark | Stereon | OCCT | Ratio | Spread | Target | Met |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ADR-0004 | Compact() after 50 deletions, 10⁴ faces | 6.4 ms | — | — | 4% | report | — |
| ADR-0004 | Compact() after 50 deletions, 10⁶ faces | 1101 ms | — | — | 7% | report | — |
| ADR-0005 | 100 edits without Compact(), 10⁴ faces | 281 ms | — | — | 15% | report | — |
| ADR-0005 | 100 edits without Compact(), 10⁶ faces | 48049 ms | — | — | 3% | report | — |
| B0.3 | OCCT's ancestor map alone, 10⁴ faces | — | 7.22 ms | — | 10% | report | — |
| B0.3 | OCCT's ancestor map alone, 10⁶ faces | — | 2648 ms | — | 9% | report | — |
| B0.3 | Edges of every vertex, 10⁴ faces | 502 µs | 13.5 ms | 27 | 23% | report | — |
| B0.3 | Edges of every vertex, 10⁶ faces | 53 ms | 3158 ms | 60 | 10% | ≥ 5× | yes |
| B0.3 | Faces of every edge, 10⁴ faces | 956 µs | 7.8 ms | 8.2 | 8% | report | — |
| B0.3 | Faces of every edge, 10⁶ faces | 123 ms | 2516 ms | 20 | 14% | ≥ 5× | yes |
| B0.3 | Build the sheet, 10⁴ faces | 9.53 ms | 231 ms | 24 | 3% | report | — |
| B0.3 | Build the sheet, 10⁶ faces | 1222 ms | 23045 ms | 19 | 39% | report | — |
| B0.5 | Delete one face, 10⁴ faces | 2.8 ms | 15.9 ms | 5.7 | 4% | report | — |
| B0.5 | Delete one face, 10⁶ faces | 475 ms | 1961 ms | 4.1 | 8% | report | — |
| B0.5 | Replace one face's surface, 10⁴ faces | 2.8 ms | 16.3 ms | 5.8 | 3% | report | — |
| B0.5 | Replace one face's surface, 10⁶ faces | 486 ms | 1962 ms | 4.0 | 5% | report | — |
| B1.1 | Curve D0, cubic, 24 poles | 51.4 ns | 84.8 ns | 1.65 | 2% | ≥ parity | yes |
| B1.1 | Curve D0, rational cubic | 57.3 ns | 88.7 ns | 1.55 | 0% | ≥ parity | yes |
| B1.1 | Curve D0, degree 7 | 94.7 ns | 152 ns | 1.61 | 1% | ≥ parity | yes |
| B1.1 | Curve D0–D2, cubic | 94.6 ns | 98.6 ns | 1.04 | 1% | ≥ parity | parity |
| B1.1 | Curve D0–D2, rational cubic | 104 ns | 146 ns | 1.40 | 1% | ≥ parity | yes |
| B1.1 | Surface D0, bicubic 12 × 12 | 109 ns | 173 ns | 1.59 | 1% | ≥ parity | yes |
| B1.1 | Surface D0, rational bicubic | 115 ns | 198 ns | 1.73 | 1% | ≥ parity | yes |
| B1.1 | Surface D0–D1, bicubic | 192 ns | 254 ns | 1.33 | 3% | ≥ parity | yes |
| B1.3 | Curve projection, rational cubic | 5.49 µs | 21 µs | 3.8 | 3% | ≥ parity | yes |
| B1.3 | Surface projection, bicubic | 103 µs | 129 µs | 1.25 | 2% | ≥ parity | yes |
| B1.4 | Approximation, 1000 points, 1e-3 mm | 1.94 ms | 35.9 ms | 18 | 17% | ≥ parity at equal tolerance | see note (not like for like) |
| B1.4 | Interpolation, 1000 points, cubic | 84.1 µs | 111 µs | 1.32 | 22% | ≥ parity | yes |
| B1.5 | Bezier decomposition, 21 segments | 4.9 µs | 7.07 µs | 1.44 | 9% | ≥ parity | yes |
| B1.5 | 100 knot insertions | 87.8 µs | 228 µs | 2.6 | 1% | ≥ parity | yes |
| B1.6 | Arc length of the whole curve, rational cubic | 106 µs | 113 µs | 1.07 | 1% | ≥ parity at equal accuracy | parity |
| B1.6 | Frenet frame and curvature, rational cubic | 135 ns | 204 ns | 1.51 | 1% | ≥ parity | yes |
| B1.6 | Parameter at a given length | 182 µs | 293 µs | 1.61 | 1% | ≥ parity | yes |
| B1.6 | Principal curvatures and directions, bicubic | 450 ns | 405 ns | 0.90 | 2% | ≥ parity | no |

## Memory and copying

Median over the repetitions of each counter. Live heap is the bytes the sheet holds in live allocations with the shape alive, measured the same way for both libraries (bench/Memory.hpp); pages copied are the pages of the edited snapshot not shared with the body (ADR-0004, ADR-0005). Stereon / OCCT below 100 % means Stereon uses less.

| ID | Measure | Stereon | OCCT | Stereon / OCCT | Target | Met |
| --- | --- | --- | --- | --- | --- | --- |
| B0.4 | Live heap per face, 10⁴ faces | 1,270 B | 1,611 B | 79% | report | — |
| B0.4 | Live heap per edge, 10⁴ faces | 629 B | 797 B | 79% | report | — |
| B0.4 | Table pages per face, 10⁴ faces | 1,267 B | — | — | report | — |
| B0.5 | Pages copied, one surface replaced, 10⁴ faces | 0.337% | — | — | report | — |
| B0.5 | Heap added, one surface replaced, 10⁴ faces | 52,736 B | 479,952 B | 11% | report | — |
| B0.5 | Pages copied, one face deleted, 10⁴ faces | 1.347% | — | — | report | — |
| B0.5 | Heap added, one face deleted, 10⁴ faces | 104,128 B | 479,856 B | 22% | report | — |
| ADR-0005 | Body and its edited copy over the body, 100 edits, 10⁴ faces | 1.451× | — | — | report | — |
| B0.4 | Live heap per face, 10⁶ faces | 1,255 B | 1,601 B | 78% | < 50% of OCCT | no |
| B0.4 | Live heap per edge, 10⁶ faces | 627 B | 800 B | 78% | report | — |
| B0.4 | Table pages per face, 10⁶ faces | 1,253 B | — | — | report | — |
| B0.5 | Pages copied, one surface replaced, 10⁶ faces | 0.003% | — | — | < 1% of body | yes |
| B0.5 | Heap added, one surface replaced, 10⁶ faces | 983,168 B | 48,000,432 B | 2% | report | — |
| B0.5 | Pages copied, one face deleted, 10⁶ faces | 0.014% | — | — | < 1% of body | yes |
| B0.5 | Heap added, one face deleted, 10⁶ faces | 1,034,560 B | 48,000,048 B | 2% | report | — |
| ADR-0005 | Body and its edited copy over the body, 100 edits, 10⁶ faces | 1.008× | — | — | < 1.2× | yes |

## Notes

- **B1.6/ArcLength/rational-cubic**, `length_mm`: 111.80161013322777 (Stereon), 111.80161013322775 (OCCT).
- **Approximation is not like for like.** Stereon returns the fewest-pole least-squares cubic within the budget; OCCT's `GeomAPI_PointsToBSpline` also optimises parameters and continuity. Read the ratio as indicative.
