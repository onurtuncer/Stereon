# Phase 1 benchmark report — curves and surfaces (B1.1–B1.6)

## Setup

| | |
| --- | --- |
| Date | 2026-10-19 |
| Commit | 9c0be02be1adbb75aeb7b90932cde15e5e6f6281 |
| Ref | main |
| Workflow run | https://github.com/onurtuncer/Stereon/actions/runs/39236132501 |
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
| Repetitions | 15 |

## Results

Median CPU time over the repetitions, both libraries interleaved: per evaluated parameter, projected point, property or length query; per fit (B1.4), per 100 insertions and per decomposition (B1.5). Ratio = OCCT / Stereon, so > 1 means Stereon is faster; within 10 % of 1 is parity. Spread is the largest relative deviation of any repetition from its median, over both libraries.

| ID | Benchmark | Stereon | OCCT | Ratio | Spread | Target | Met |
| --- | --- | --- | --- | --- | --- | --- | --- |
| B1.1 | Curve D0, cubic, 24 poles | 51.4 ns | 84.9 ns | 1.65 | 4% | ≥ parity | yes |
| B1.1 | Curve D0, rational cubic | 57.3 ns | 88.3 ns | 1.54 | 3% | ≥ parity | yes |
| B1.1 | Curve D0, degree 7 | 94.6 ns | 153 ns | 1.61 | 3% | ≥ parity | yes |
| B1.1 | Curve D0–D2, cubic | 94.7 ns | 98.7 ns | 1.04 | 3% | ≥ parity | parity |
| B1.1 | Curve D0–D2, rational cubic | 104 ns | 139 ns | 1.33 | 1% | ≥ parity | yes |
| B1.1 | Surface D0, bicubic 12 × 12 | 109 ns | 172 ns | 1.58 | 1% | ≥ parity | yes |
| B1.1 | Surface D0, rational bicubic | 115 ns | 198 ns | 1.73 | 1% | ≥ parity | yes |
| B1.1 | Surface D0–D1, bicubic | 191 ns | 255 ns | 1.33 | 2% | ≥ parity | yes |
| B1.3 | Curve projection, rational cubic | 5.5 µs | 20.5 µs | 3.7 | 3% | ≥ parity | yes |
| B1.3 | Surface projection, bicubic | 104 µs | 128 µs | 1.24 | 3% | ≥ parity | yes |
| B1.4 | Approximation, 1000 points, 1e-3 mm | 1.95 ms | 36.7 ms | 19 | 31% | ≥ parity at equal tolerance | see note (not like for like) |
| B1.4 | Interpolation, 1000 points, cubic | 84.1 µs | 113 µs | 1.34 | 16% | ≥ parity | yes |
| B1.5 | Bezier decomposition, 21 segments | 4.97 µs | 7.3 µs | 1.47 | 10% | ≥ parity | yes |
| B1.5 | 100 knot insertions | 88 µs | 227 µs | 2.6 | 1% | ≥ parity | yes |
| B1.6 | Arc length of the whole curve, rational cubic | 106 µs | 112 µs | 1.06 | 3% | ≥ parity at equal accuracy | parity |
| B1.6 | Frenet frame and curvature, rational cubic | 135 ns | 189 ns | 1.41 | 2% | ≥ parity | yes |
| B1.6 | Parameter at a given length | 182 µs | 293 µs | 1.61 | 2% | ≥ parity | yes |
| B1.6 | Principal curvatures and directions, bicubic | 442 ns | 409 ns | 0.93 | 3% | ≥ parity | parity |

## Notes

- **B1.6/ArcLength/rational-cubic**, `length_mm`: 111.80161013322777 (Stereon), 111.80161013322775 (OCCT).
- **Approximation is not like for like.** Stereon returns the fewest-pole least-squares cubic within the budget; OCCT's `GeomAPI_PointsToBSpline` also optimises parameters and continuity. Read the ratio as indicative.
