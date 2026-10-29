# Benchmark reference run — topology (B0.3–B0.5), curves and surfaces (B1.1–B1.6)

## Setup

| | |
| --- | --- |
| Date | 2026-10-20 |
| Commit | 0e3b9848dadf62ccbaee57c9ac41da4d0244f510 |
| Ref | topo/typed-pools |
| Workflow run | https://github.com/onurtuncer/Stereon/actions/runs/39437110770 |
| Runner | GitHub Actions ubuntu-24.04 (Linux, X64), container gcc:16.2 |
| CPU | INTEL(R) XEON(R) PLATINUM 8573C |
| CPUs | 4 |
| Memory | 15.6 GB |
| Kernel | Linux 6.17.0-1022-azure |
| Compiler | g++ (GCC) 16.2.0 |
| OCCT tag | V8_0_1 |
| OCCT commit | b8f597c677811d1f9f4d8a97f5ae2825c0353a42 |
| OCCT compiler | c++ (GCC) 16.2.0 |
| OCCT build type | Release |
| Repetitions | 3 |
| memory_measure | live heap, glibc mallinfo2 (bench/Memory.hpp) |

## Results

Median CPU time over the repetitions, both libraries interleaved: per evaluated parameter, projected point, property or length query; per fit (B1.4), per 100 insertions and per decomposition (B1.5); per build, per query over the whole sheet or per edit (B0.3, B0.5, ADR-0004, ADR-0005). Ratio = OCCT / Stereon, so > 1 means Stereon is faster; within 10 % of 1 is parity, and a row with a ratio target of its own (B0.3: ≥ 5×) is met when it reaches it. A target of "report" is not a gate. Spread is the largest relative deviation of any repetition from its median, over both libraries.

| ID | Benchmark | Stereon | OCCT | Ratio | Spread | Target | Met |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ADR-0004 | Compact() after 50 deletions, 10⁴ faces | 5.33 ms | — | — | 1% | report | — |
| ADR-0004 | Compact() after 50 deletions, 10⁶ faces | 819 ms | — | — | 21% | report | — |
| ADR-0005 | 100 edits without Compact(), 10⁴ faces | 273 ms | — | — | 13% | report | — |
| ADR-0005 | 100 edits without Compact(), 10⁶ faces | 44549 ms | — | — | 3% | report | — |
| B0.5 | Delete one face, 10⁴ faces | 2.74 ms | 29 ms | 11 | 9% | report | — |
| B0.5 | Delete one face, 10⁶ faces | 445 ms | 2969 ms | 6.7 | 2% | report | — |
| B0.5 | Replace one face's surface, 10⁴ faces | 2.74 ms | 27.7 ms | 10 | 2% | report | — |
| B0.5 | Replace one face's surface, 10⁶ faces | 450 ms | 2984 ms | 6.6 | 40% | report | — |
| B1.1 | Curve D0, cubic, 24 poles | 24.6 ns | 68.9 ns | 2.8 | 1% | ≥ parity | yes |
| B1.1 | Curve D0, rational cubic | 31.1 ns | 74.3 ns | 2.4 | 2% | ≥ parity | yes |
| B1.1 | Curve D0, degree 7 | 59.8 ns | 128 ns | 2.1 | 1% | ≥ parity | yes |
| B1.1 | Curve D0–D2, cubic | 65.6 ns | 76.2 ns | 1.16 | 0% | ≥ parity | yes |
| B1.1 | Curve D0–D2, rational cubic | 71.1 ns | 114 ns | 1.60 | 1% | ≥ parity | yes |
| B1.1 | Surface D0, bicubic 12 × 12 | 69 ns | 132 ns | 1.91 | 2% | ≥ parity | yes |
| B1.1 | Surface D0, rational bicubic | 77.1 ns | 161 ns | 2.1 | 1% | ≥ parity | yes |
| B1.1 | Surface D0–D1, bicubic | 147 ns | 207 ns | 1.41 | 5% | ≥ parity | yes |

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
| ADR-0005 | Body and its edited copy over the body, 100 edits, 10⁴ faces | 1.324× | — | — | report | — |
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
| ADR-0005 | Body and its edited copy over the body, 100 edits, 10⁶ faces | 1.006× | — | — | < 1.2× | yes |
