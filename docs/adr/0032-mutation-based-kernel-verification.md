# ADR-0032: Mutation-Based Verification of Kernel Contracts

- **Status:** Accepted
- **Date:** 2026-10-16
- **Revised:** 2026-10-17 (scheduled campaigns run weekly)
- **Accepted:** 2026-10-17 by Onur Tuncer, on the Phase 0 pilot evidence in [docs/planning/mutation-pilot.md](../planning/mutation-pilot.md)
- **Owner:** Onur Tuncer
- **Phase:** 0 (verification policy); applied incrementally as kernel capabilities mature
- **Related:** ADR-0002, ADR-0004, ADR-0006, ADR-0007, ADR-0016, ADR-0018, ADR-0026–ADR-0031

## Context

Line and branch coverage demonstrate execution, not whether tests can detect incorrect behavior. A CAD kernel can execute every line of an algorithm while silently misclassifying an intersection, breaking an incidence relation, accepting a stale handle, or losing an operation's entity lineage. In particular, a validity checker may approve a broken result if both producer and checker share the same mistaken assumption.

Stereon already has Catch2 unit tests, RapidCheck property tests, contract probes, differential tests against OCCT, and coverage instrumentation. Those measures are complementary, but none measures how sensitive the tests are to plausible implementation faults. ADR-0018 requires a seeded-fault corpus for the independent validity checker. We need a consistent strategy for testing the tests across the kernel, including observable contracts that ordinary source-level mutation operators do not model.

Here, **mutation testing** means deliberately introducing a temporary fault into *implementation code* and asking whether the unmodified tests detect it. It is distinct from *topological mutation*, the legitimate edit operations specified by ADR-0027. Mutants must never be shipped or committed as production behavior.

## Options

1. **Rely on coverage, property tests, fuzzing and differential tests alone.** Lower cost, but provides no systematic evidence that the assertions distinguish correct from faulty behavior.
2. **Require one mutation-testing framework and a universal mutation-score threshold.** Easy to describe, but brittle across C++26 modules, compilers, equivalent mutations, and differing fault severity.
3. **Use a tool-independent, risk-based mutation-verification policy with conventional and CAD-specific fault operators.** Requires a catalogue, triage and dedicated CI capacity, but connects failure modes directly to architectural guarantees.

## Decision

**Adopt option 3. Stereon shall evaluate the fault-detection effectiveness of its tests through controlled mutation experiments.** The verification strategy shall combine ordinary source-level mutants with curated, domain-specific faults derived from the kernel's ADR contracts. A high code-coverage value alone shall not establish adequacy of tests for a supported kernel operation.

This ADR sets a verification policy, not a dependency on a particular mutation framework, compiler, or score target. Its requirements are phased with the corresponding kernel functionality.

### 1. Two complementary mutation families

**Conventional source-level mutants** include relational and arithmetic operator substitutions, condition negation, Boolean-return changes, and removal of selected checks or state updates. Only mutations that preserve a meaningful build and exercise an observable contract are scoreable; a compile failure is not evidence of test adequacy.

**Domain-specific mutants** model plausible CAD-kernel defects and must be explicitly linked to the behavioral obligation being challenged:

| Domain and governing ADR | Example seeded fault | Required independent observation |
| --- | --- | --- |
| Robust predicates — ADR-0002 | Reverse an orientation sign or mishandle a filtered-predicate fallback | An exact/independent geometric oracle rejects the classification |
| Tolerance policy — ADR-0001/0028 | Use an angular or solver tolerance as a length tolerance; change a boundary comparison | Dimensional/boundary tests reject the decision |
| Handle identity — ADR-0004 | Accept a stale generation or a foreign-store handle | Checked-handle or API contract probe detects misuse |
| Topology / Euler operators — ADR-0027 | Skip the reciprocal coedge update or leave a radial ring open | Independent structural traversal reports the violated rule |
| Validity checker — ADR-0018 | Suppress a specified invalidity finding | A known-invalid seeded shape is no longer accepted by the test suite |
| Operations / failure semantics — ADR-0006/0026 | Publish a partial result after an injected failure | Rollback, immutable-input and result-contract tests detect it |
| Intersections — ADR-0029 | Report a tangency as transverse or coerce an unresolved result to disjoint | Geometric oracle and outcome-specific assertions detect it |
| Lineage — ADR-0030 | Omit one split descendant or mark a generated entity as a continuation | Manifest-completeness and ownership checks detect it |
| Attributes — ADR-0031 | Silently copy a unique ID across a split or discard a merge conflict | Policy-specific propagation tests detect it |
| Determinism — ADR-0007 | Make a result depend on traversal/thread ordering | G1 comparison across permitted thread counts detects it |

A mutant shall identify one conceptual fault and the affected contract. Complex coordinated faults may be useful later but shall not obscure which requirement a test proves.

### 2. Oracle independence and test meaning

- The **unmutated** implementation and its relevant test suite must pass before mutant results are considered.
- Mutate production code, not the test assertions or reference oracle used to judge that mutation. The test runner and reporting mechanism must remain trustworthy.
- When testing an operation that uses a validity checker, do not rely solely on that same checker's approval: use independently derived invariants, analytic examples, differential reference results, or a separately implemented checker traversal as appropriate to the guarantee.
- When mutating the validity checker itself, use a curated corpus containing **known-invalid** and **known-valid** shapes with specified expected findings. A mutant that weakens rejection cannot be killed merely by another call through the same defective rule.
- Model numerical uncertainty honestly: tolerances, valid domains, exact-predicate requirements and explicit unresolved outcomes remain governed by ADR-0001/0002/0028/0029. Test assertions shall not impose an unjustified exact answer for genuinely ambiguous geometry.
- External invalid data must still be handled as an error under ADR-0006/0016; do not turn a programmer-contract termination into an acceptable substitute for a public API's structured error.

### 3. Classification, evidence and metric

Each experiment records the source revision, mutant ID, operator, location, affected ADR/contract, test selection, seed/corpus input where relevant, compiler/build preset, result and diagnostic evidence. Preserve reproducible surviving cases for regression-test development.

Use the following mutually distinguishable outcomes:

- **Killed:** the previously passing test suite reproducibly fails because of the mutation, with diagnostic evidence attributable to it.
- **Survived:** an executable, exercised and non-equivalent mutant passes the relevant tests; a verification gap for triage.
- **Uncovered:** an executable mutant is not reached by the selected tests; report separately and treat as an adequacy gap, not a kill.
- **Equivalent:** expert review shows no observable difference within the declared supported semantics; retain the rationale and exclude from the denominator.
- **Invalid:** the mutant cannot be built or executed due to mutation mechanics rather than the property under test; exclude from the denominator.
- **Inconclusive:** infrastructure failure, non-reproducible timeout or flaky outcome; rerun or triage, never count as killed.

The **mutation score** is `killed / (killed + survived + uncovered)` over valid, classified, non-equivalent mutants in the stated scope. Always publish numerator, denominator, exclusions, operator set and scope with any score; do not compare scores drawn from incompatible mutation sets as though they were equivalent. Mutant sampling, if used, must be explicitly disclosed.

A surviving mutant is not automatically a kernel defect: it is a test-effectiveness finding requiring an additional test, a revised contract, a justified equivalence decision, or an explicit waiver. Waivers require owner, reason, affected contract and periodic review.

### 4. Toolchain and ownership boundary

ADR-0016 remains authoritative: **GCC 16+ is the reference compiler and C++26/module/contract requirements are not relaxed for a mutation tool**.

An LLVM/Clang-based framework (for example Mull) may be piloted on components that it can build faithfully. It is not a mandatory dependency. If a tool cannot reproduce Stereon's reference semantics, use a controlled source-level mutation harness compiled and executed with the reference GCC preset for critical contract faults. Experimental results from a reduced-feature Clang build shall be labelled as such and shall not substitute for reference-build verification.

Mutations run in disposable worktrees/build directories. Tool-generated mutants, temporary altered modules and instrumented binaries must not be merged or published as normal build artefacts. Keep the original tests, fixtures and contract configuration fixed across each valid baseline/mutant comparison; make any unavoidable differences explicit.

### 5. CI and resource policy

1. **Phase 0 pilot:** begin with bounded campaigns in `core/`, `robust/` and simple `curve/` algorithms; measure tool compatibility and the cost of rebuilding C++ modules. Include the stale-handle, orientation and numerical-boundary cases.
2. **Incremental development:** run selective mutations for changed, compilable components where instrumentation is reliable. Initially produce non-blocking reports rather than imposing an arbitrary score threshold.
3. **Scheduled/full campaigns:** run broader operator and domain-fault sets off the critical PR path, retaining the mutant inventory and reproducible survivor evidence. They run **weekly**, and on request, on the reference GCC 16 build (the `Mutation` workflow); a week in which the campaign does not run, or its baseline fails, is itself a finding to triage.
4. **Release gate:** once a domain-fault suite is validated for an implemented capability, every designated **critical seeded fault** for that capability must be detected by its specified independent oracle. Ordinary mutation scores are reviewed by subsystem and trend, not held to a universal 100% threshold.
5. **Budget control:** use explicit per-campaign time limits, targeted tests and deterministic sampling. Report unexecuted mutants; do not silently classify a budget-exhausted run as success.

No additional CI job, tool integration, or full mutation result is claimed by this ADR. Those are implementation tasks subject to feasibility and measured CI cost.

## Consequences

- Test reviews gain evidence that assertions actually distinguish incorrect geometry and topology, rather than merely execute code.
- Architecture decisions can be traced through a contract → fault model → test/oracle → detection record. Failures can be prioritized by severity and contract criticality.
- Mutation experiments can expose common-mode errors in validity checks and independent reference implementations.
- Compiler/module compatibility, equivalent-mutant triage and corpus maintenance add expense. Mutation scores can also be gamed; the curated critical-fault inventory and oracle independence are more important than one aggregate number.
- Mutation testing complements, and does not replace, property-based tests, fuzzing, differential testing, sanitizers, contracts or coverage.

## Verification

Before this ADR advances from Proposed to Accepted, demonstrate a documented Phase 0 pilot:

1. An unmodified reference GCC 16 test baseline passes; a repeatable source-level mutation harness generates executable mutants without changing test/oracle source or the production branch.
2. Run a curated set covering at least one boundary comparison, one orientation/predicate fault, one stale-handle acceptance fault and one independent-checker fault (the last may use the ADR-0018 seeded-fault harness). Each implemented, designated **critical** fault is detected with the expected diagnostic.
3. Intentionally add an ordinary mutant that survives an initially inadequate test; add an assertion that kills it, and show the difference between code coverage and mutation effectiveness.
4. Show correct reporting for a killed, survived or uncovered, invalid and equivalent/inconclusive example, including the explicitly defined denominator and exclusion rationale.
5. Record the build/compiler combination, CI time and C++ module compatibility limits; confirm the reference GCC test suite is unaffected and mutants are not present in published artefacts.

Subsequent capability gates extend the curated fault catalogue to topology operations, intersections, lineage and attributes as those features become supported. A failed critical seeded fault blocks that capability's verification gate **after** its corresponding suite has been validated and adopted.

## Relationships

- [ADR-0002](0002-numeric-robustness.md): robust predicates and assurance.
- [ADR-0004](0004-topology-storage-and-handles.md): checked handles and topology.
- [ADR-0006](0006-error-handling.md): error and failure semantics.
- [ADR-0007](0007-threading-and-determinism.md): determinism.
- [ADR-0016](0016-cpp26-adoption-and-contracts.md): reference toolchain and programmer contracts.
- [ADR-0018](0018-validity-checking.md): independent validity checks and seeded-fault corpus.
- [ADR-0026](0026-modelling-operations-by-contracts.md): behavioral modelling-operation contracts.
- [ADR-0027](0027-controlled-topology-mutation-and-euler-operators.md): controlled *model* mutation, distinct from *test* mutation.
- [ADR-0028](0028-purpose-specific-numerical-tolerances.md): numerical policy.
- [ADR-0029](0029-intersection-and-classification-outcomes.md): explicit outcome classes.
- [ADR-0030](0030-entity-lineage-in-operation-results.md): lineage completeness.
- [ADR-0031](0031-attribute-propagation-and-conflicts.md): attribute semantics.

This proposal supplements these ADRs and supersedes none.
