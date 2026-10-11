# ADR-0022: Function complexity limits

- **Status:** Proposed
- **Date:** 2026-10-11
- **Owner:** Onur Tuncer
- **Phase:** 0

## Context

ADR-0015 rests on people reading code: Tier C changes (`intersect`, `boolean`, `blend`) get line-by-line review, and a Tier B reviewer must be able to explain the algorithm. Both fail once a function grows past what a reader can hold in mind. Geometric algorithms grow exactly that way: surface–surface intersection, Boolean classification and blending accumulate degenerate-case branches one bug report at a time, and a function that handles twelve cases in one body passes the obvious tests while hiding the failures on hard cases that ADR-0015 warns about.

Review judgement alone does not hold a limit, because each addition is small. A numeric limit, checked in CI, does. The limit also shapes the code towards how the project verifies it: one function per geometric case is what the differential tests, the corpus and the history records of ADR-0011 can point at.

Tooling constrains the choice. clang-tidy cannot parse the GCC 16 module and reflection build (see `.clang-tidy`), so its `readability-function-cognitive-complexity` check does not run; the CI lint job uses Python checks that read the sources instead (`tools/check_layering.py`, `tools/check_epsilon.py`). `lizard` measures cyclomatic complexity and function length the same way, without compiling, and handles module declarations, contract specifiers and expansion statements in this codebase. At the time of writing the tree has 302 functions with an average cyclomatic complexity of 1.4 and a maximum of 5.

## Options

1. **No numeric limit; review judgement only.** Costs nothing; erodes one branch at a time.
2. **Cyclomatic complexity per function (McCabe), checked with `lizard` in the lint job.** Counts decision points; cheap to measure on this build; blind to nesting, so a flat `switch` and a four-deep loop nest with the same count read very differently.
3. **Cognitive complexity per function (SonarSource's definition), checked with clang-tidy.** Penalises nesting, which is what makes geometric case analysis hard to read; needs a clang-tidy that can parse the build, which none does today.
4. **Length limits only.** Catches sprawl, not branching; a 40-line function can still have twenty paths.

## Decision

Adopt option 2 together with a length limit (option 4) now, and add option 3 when a clang-tidy can parse the build. Limits follow the ADR-0015 tiers:

| Scope | Cyclomatic complexity | Non-comment lines | Cognitive complexity (when measurable) |
| --- | --- | --- | --- |
| Tier C: `intersect`, `boolean`, `blend` | ≤ 10 | ≤ 60 | ≤ 15 |
| Tier B: `robust`, `topo`, `build`, `mesh` | ≤ 12 | ≤ 80 | ≤ 20 |
| Tier A: `core`, `curve`, `surface`, `io`, `python`, `bench`, the viewer | ≤ 15 | ≤ 100 | ≤ 25 |
| `tests/` | ≤ 15 | ≤ 150 | ≤ 25 |

- **What is measured.** Every function, member function and function template, with the bodies of lambdas defined inside it counted as part of it. Moving a block into a lambda called once therefore changes nothing; splitting it into a named function does.
- **How to split.** A function over the limit is split along its case analysis, not along line numbers: one function per contact type of ADR-0010, per Boolean classification outcome of ADR-0013, per file-format chunk of ADR-0012. The names record the case (ADR-0017), so the decomposition documents the algorithm.
- **Exemptions are explicit.** A function that legitimately exceeds a limit carries the comment `// complexity-exempt: <reason>` on the line before it. Expected reasons: a `switch` over a closed enum where every case is one line, a parser state machine in `io`, or a reference implementation kept deliberately literal for comparison. The exemption is reviewed like code, and the check reports the number of exemptions in every run.
- **Generated code** needs no exemption: tables built with expansion statements and reflection (ADR-0016) contribute few branches in source, which is one more reason to generate dispatch rather than write it.
- **Enforcement.** `tools/check_complexity.py` wraps `lizard`: it lists the tracked `.cppm`, `.cpp` and `.hpp` files itself (lizard skips `.cppm` otherwise), strips `STEREON_PRE`, `STEREON_POST` and their `HOT` forms before measuring so that functions report under their own names, applies the table above by directory, honours the exemption marker, and fails on any violation. The lint job runs it with a pinned lizard version, as it pins clang-format, and CTest runs it like the other tools. Thresholds live in the script next to the directory map, in one place.
- **Tests are not exempt.** A test case is a list of assertions and so is allowed to be long, but a test with branching logic of its own is a test that can be wrong; the complexity limit holds there too.

## Consequences

- Tier C algorithms are forced into case-by-case decomposition, which is also how reviewers, differential tests and corpus failures reason about them.
- Some algorithms meet the limit naturally only with effort: a Newton iteration with several exit conditions, a tessellator's edge-splitting loop, a STEP parser. They are decomposed or carry a reviewed exemption; neither happens silently.
- Cyclomatic complexity is blind to nesting until the cognitive measure can run; reviewers watch nesting depth in the meantime.
- The lint job gains a pinned Python dependency.
- The number of small functions grows. Accepted: with ADR-0017 names they read as the algorithm's outline.
- A limit that is reviewed away through exemptions is no limit; the exemption count is therefore reported, and a rising count is a review item at each phase gate.

## Verification

- `tools/check_complexity.py` exits non-zero on any violation; the lint job and CTest run it. Pass at adoption: zero violations and zero exemptions over the current tree (302 functions, maximum complexity 5).
- The script has a `--self-test`, like `check_epsilon.py`: a synthetic function of complexity 11 under `src/boolean/` must be flagged, the same function under `src/core/` must pass, and a flagged function with the exemption marker must pass and be counted.
- Each run prints the number of exemptions; the Phase 0 exit review records it.
- When a clang-tidy job exists, `readability-function-cognitive-complexity` enforces the last column, and this ADR is revised with the first measurements.
