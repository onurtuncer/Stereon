# ADR-0032 Phase 0 mutation pilot

**Date:** 2026-10-17
**Revision measured:** `f06b521` (branch `verify/mutation-pilot`)
**Status:** evidence for ADR-0032's Verification section, on which the ADR was accepted on 2026-10-17 with weekly campaigns.

This report records the pilot that ADR-0032 requires before it can move from Proposed to Accepted. Each section answers one item of the ADR's Verification list. The harness is `tools/mutate.py`, the catalogue is `tools/mutation/pilot.json`, and the method is documented in `docs/testing.rst` (*Mutation testing*).

## 1. Baseline and harness

- Reference toolchain: g++ 16.2.0 (MSYS2 UCRT64), CMake 4.4.4, preset `gcc16-debug` (contracts enforced, checked handles on), Windows 11. The `Mutation` workflow runs the same harness in the `gcc:16.2` Linux container.
- Each run checks out the commit into a disposable git worktree under `build/mutation`, configures it, and builds and runs every selected test set unmodified. The baseline passed before any mutant was applied.
- Each mutant is one unique source replacement, applied alone. The harness rebuilds incrementally, runs the selected tests and restores the file with `git checkout`. Tests and oracles are never mutated, and the worktree is removed at the end. The developer's tree and branches are never touched.

## 2. Critical faults

All five designated critical faults were detected, each with its declared diagnostic:

| Mutant | Seeded fault | Contract | Detected by |
| --- | --- | --- | --- |
| PRED-FALLBACK | The filtered `Orient2D` trusts the floating-point sign and skips the exact fallback | ADR-0002 | `robust: adaptive and exact stages agree on nearly degenerate floating-point input`, `predicates are exact on tiny and huge coordinates`, `predicates survive the adversarial sweep` |
| PRED-SIGN | `Orient2D` returns the opposite sign | ADR-0002 | compile-time tests: `static_assert`s in `PredicateTests.cpp` reject the build |
| HANDLE-STALE | A handle's generation is not compared: stale handles are accepted | ADR-0004 | `core: generation tables detect stale handles`, `generation tables keep their invariants under random use` |
| BOUND-INTERVAL | `Contains(interval, x)` treats the lower bound as open | ADR-0001, ADR-0028 | compile-time tests in `BoxTests.cpp` and `ContextTests.cpp` |
| CHECK-KNOTS | The knot-vector check stops rejecting an interior knot repeated degree + 1 times | ADR-0016, standing in for ADR-0018 | `curve: MakeBSplineCurve validates data from outside the kernel` aborts when the expected error is missing |

**About the independent-checker fault.** ADR-0018's validity checker does not exist yet; topology is still a stub. CHECK-KNOTS uses the nearest independent checker the kernel has: the checked factories' validation, which every file read (`ReadGeometry`) and every outside input passes through. It is judged against the factory tests' corpus of known-invalid inputs. Once ADR-0018's checker and its seeded-fault corpus exist, they should replace this case.

## 3. Coverage is not mutation effectiveness

An earlier campaign classified BOUND-FRAME (in `WithinFrameAllowance`, `std::abs(x) <= FrameAllowance` changed to `<`) as **survived**. The harness had planted an execution trap at the same site, and the curve tests failed on it, so they do reach the check: line coverage there is complete.

The only test near the boundary used an error of 1e-14, about a hundred times inside the allowance, so nothing pinned the boundary itself.

`tests/curve/EllipseTests.cpp` now has "the frame allowance is a closed bound". A frame whose error is exactly `FrameAllowance` (2⁻⁴⁰, exact in binary) must be accepted, and the next double above it must be rejected. In the final campaign the same mutant is **killed** by that test.

## 4. Reporting every outcome

The final campaign, on revision `f06b521`:

| Outcome | Count | Mutants |
| --- | --- | --- |
| Killed | 8 | the five critical faults; BOUND-FRAME; ARITH-BASIS (`+` to `-` in the Cox–de Boor recurrence); REL-SPAN (`>=` to `>` in the span search) |
| Survived | 0 | (BOUND-FRAME survived in the first campaign; see section 3) |
| Uncovered | 0 | (in the first campaign BOUND-FRAME was uncovered under a selection that left out the tests reaching it; the catalogue now selects the curve tests) |
| Equivalent | 1 | EQUIV-PRUNE (`>` to `>=` in best-first pruning), with an expert rationale in the catalogue |
| Invalid | 1 | INVALID-DEMO: production code that does not build |
| Inconclusive | 1 | HANG-DEMO: the span search never terminates, and the test timed out |

**Score:** killed / (killed + survived + uncovered) = 8 / 8 = 1.00 over this catalogue.
- Excluded from the denominator: EQUIV-PRUNE (equivalent, with rationale), INVALID-DEMO (invalid) and HANG-DEMO (inconclusive).
- The operator set and scope are exactly the eleven mutants of `tools/mutation/pilot.json`, so this score says nothing about code outside them. ADR-0032 forbids comparing it with scores from other mutation sets.

Two classification rules came out of the first runs and are now part of the harness:
- **A `static_assert` failing in a test file is a kill** by a compile-time test. Only a build failure in production code is Invalid. PRED-SIGN and BOUND-INTERVAL were first misreported as Invalid.
- **A test selection that matches no test is Inconclusive,** never a pass, and the baseline refuses such a selection.

## 5. Cost and limits

- **Time:** about 30 minutes for the eleven mutants on a 4-core laptop. The report's 41,911 s wall time includes an overnight suspend during HANDLE-STALE's rebuild: its log shows 5 s of tests and no output between 23:44 and 10:55. Per mutant, 20 s (curve changes) to 6 min. The cost is the rebuild: a change in `core` rebuilds every module that imports it, so core mutants cost the most, and C++ modules make the rebuild order strict.
- **Tools:** no LLVM-based framework (such as Mull) was used. It cannot build the reference configuration: GCC 16 modules, contracts and reflection (ADR-0016). The harness is a source-level one on the reference compiler, as ADR-0032 section 4 allows.
- **Isolation:** nothing from a mutant reaches a branch or an artefact. The worktree is detached and removed, and its build directory sits inside it. Saved logs (`build/mutation/logs`) and the report are the only outputs.
- **Windows:** one transient "permission denied" on a compiler intermediate (file locking by another process) needed a rerun of a local build. It did not affect a classified result.
- **CI:** the `Mutation` workflow runs weekly and on request, and fails only when the baseline fails or a critical fault goes undetected, as ADR-0032 section 5 asks. It has not run on GitHub yet: a scheduled workflow becomes active once it is on the default branch.

## Next

- Run the workflow once on GitHub after merging, and record the CI time.
- Extend the catalogue as capabilities land: topology and Euler operators (ADR-0027), the validity checker's seeded-fault corpus (ADR-0018), intersections (ADR-0029), lineage (ADR-0030) and attributes (ADR-0031).
- Add a mutant to the catalogue whenever a defect is fixed, so that its regression test is shown to catch it.
