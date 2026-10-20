# ADR-0042: Assign C++26 contracts and kernel validation distinct responsibilities

- **Project:** Stereon
- **Status:** Proposed
- **Date:** 2026-10-20
- **Owner:** Onur Tuncer
- **Phase:** 0
- **Scope:** Public input checks, internal assertions, publication gates, compile-time verification and diagnostic generation
- **Related:** Operation contracts; controlled topology mutation; topology/geometry separation; proposed executable-contract, atomic-publication and layered-validation ADRs

## Context

C++26 introduces native contract assertions through `pre`, `post` and `contract_assert`. Their evaluation semantics are implementation-defined and include three policy groups (enforcing includes both `enforce` and `quick-enforce`; record the exact semantic in build configuration):

- **Ignored:** contract predicates are not evaluated and violations do not affect execution.
- **Observing:** predicates are evaluated and violations are reported through the implementation’s observation mechanism, but execution is not required to terminate.
- **Enforcing:** violations are handled as contract failures and may terminate execution according to the implementation’s contract semantics.

Native contracts are therefore not, by themselves, a portable recoverable-input-validation mechanism. In particular, public safety and publication guarantees must not depend on observing or enforcing mode.

Stereon’s public operations must report invalid requests through their documented error policy, preserve inputs on supported failure paths, and establish required result guarantees before publication. Those promises must hold in ignored, observing and enforcing modes. A correct implementation must not rely on a native contract being evaluated, and enabling contracts must not replace or weaken ordinary validation.

Function postconditions apply on normal return. An operation that commits inside its function body cannot rely on a later native postcondition as its pre-publication gate. Native contracts also do not provide rollback, capture the old model automatically, or prove geometric correctness.

Compile-time evaluation and reflection offer complementary benefits. Small deterministic topology fixtures can be checked during compilation, and reflection can reduce repetitive diagnostic code. Neither technique infers the semantic invariants of a CAD kernel.

## Decision

**Use C++26 contracts to express and check internal programming obligations. Keep public input validation, required publication checks and atomic failure guarantees explicit and independent of contract evaluation mode.**

The following guarantees are mandatory in every mode:

- Public input validation remains active in ignored, observing and enforcing modes.
- Required pre-publication validation executes independently of native contract evaluation.
- Unsupported or unresolved input is handled through the documented recoverable error policy.
- Supported failure paths preserve the required input, aliases and externally visible state.
- Publication occurs only after all mandatory properties have been established.
- Rollback and cleanup do not depend on a contract handler, contract violation notification or stack unwinding after termination.
- Enabling observing or enforcing contracts may add diagnostics or terminate on programming defects, but must not be required to preserve public safety.
- Ignoring contracts must not permit a correct implementation to publish a result that violates a mandatory invariant.

The following responsibility allocation is normative:

|Responsibility                       |Preferred mechanism                                     |Required boundary                                                                                                              |
|-------------------------------------|--------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------|
|Invalid public input                 |Explicit checks and structured errors                   |Execute before unsafe access or mutation in every mode, using the existing error/exception policy.                             |
|Internal programming assumptions     |C++26 contracts where supported                         |Supplement checked inputs and reviewed invariants; never provide the sole enforcement of a public safety requirement.          |
|Required validation before commit    |Explicit executable predicates and valid scoped evidence|Establish every promised property before publication in every mode; violated or unresolved required properties block commit.   |
|Rollback and atomic publication      |Transactions and RAII                                   |Preserve public inputs and aliases on supported failures; cleanup must not depend on a contract handler running.               |
|Small deterministic topology fixtures|`constexpr` and `static_assert`                         |Exercise the same semantic predicates or transformations used at runtime where practical; retain runtime and adversarial tests.|
|Repetitive diagnostic machinery      |Reflection where useful                                 |Generate mechanical structure and field information; specify semantic rules and stable diagnostic identifiers explicitly.      |

### 1. Public API validation remains ordinary control flow

A public Result-returning operation or query checks handle identity, ownership, generation/liveness, parameter validity and supported input profiles before unsafe dereference or mutation. These checks execute in ignored, observing and enforcing modes. Expensive geometric requirements may be established through current, appropriately scoped validation evidence under the layered-validation policy.

Unsupported or numerically unresolved inputs follow the documented error policy. They must not be turned into contract violations merely to avoid designing recoverable errors. Public checks remain active in release builds and builds that ignore native contracts.

Illustrative code, not an existing Stereon API:

```cpp
if (!builder.Contains(edge))
    return Failure(InvalidHandle{});

if (!Supports(request.Profile))
    return Failure(UnsupportedProfile{});
```

The O(1) value-returning accessor and Context-modifier exceptions in ADR-0006 remain programmer obligations under ADR-0016. This proposal does not change those APIs or add recoverable error returns to them.

### 2. Native contracts express internal obligations

Express `pre`, `post` and `contract_assert` through the existing `STEREON_PRE`, `STEREON_POST`, `STEREON_ASSERT` and `STEREON_HOT_*` macros from ADR-0016 for internal assumptions whose violation indicates a programming defect, such as a prepared plan being applied to the expected builder revision or an internal result having its promised shape.

Contract predicates must be observationally pure with respect to kernel state and must not rely on their evaluation count. Do not place edits, rollback, cache updates required for correctness, certificate issuance, or required logging inside a predicate. Predicates should be bounded and non-throwing where practical; they must not dereference unvalidated external handles.

The explicit validator and a native assertion may share a predicate implementation. Publication correctness must not depend on the assertion executing. Independent validation remains separately required where specified by the validation ADR; predicate reuse is not evidence of validator independence.

Mode-specific behavior is limited to diagnostics and defect handling:

- In **ignored mode**, the contract contributes no runtime check. All required guarantees must already be established by ordinary control flow, explicit predicates, transactions and reviewed construction rules.
- In **observing mode**, a violation may be reported, but execution must not be assumed to stop or roll back. The operation must still be safe under the implementation’s continued execution semantics, or the violation must be treated as an internal defect whose subsequent behavior is outside the public recoverable-failure guarantee.
- In **enforcing mode**, a violation may terminate execution. It is not a recoverable public error path and must not be used as the mechanism for ordinary invalid-input handling.

An ignored contract does not authorize broken internal invariants. Correct implementations must establish their assumptions through checked inputs, valid evidence and reviewed construction rules, and must reject detected defects before publication.

### 3. Explicit postcondition checks precede commit

For an atomic operation, mandatory checks run while changes remain private or reversible in every mode:

```text
plan = validate_and_prepare(input, request)
transaction = begin_private_edit(input)
candidate = apply_edit(transaction, plan)

report = establish_required_postconditions(candidate, plan)
if report does not establish all required properties:
    return failure(report)  // transaction rolls back or discards candidate

result = prepare_complete_result(candidate)
publish(result)             // sole commit point
return success(result)
```

All fallible result preparation that could otherwise violate atomicity happens before publication. A native `post` on the enclosing function is supplementary and cannot serve as this gate after commit has occurred. Function postconditions do not check exception exits; rollback and failure guarantees require separate mechanisms and tests.

Before/after comparisons use explicitly retained counts, correspondence, snapshots or transaction metadata. This ADR does not assume a language-provided old-state capture mechanism.

If an explicit mandatory check fails, remains unresolved or reports insufficient evidence, publication is blocked regardless of contract mode. An observing notification cannot substitute for the failure path, and an enforcing violation is not the ordinary mechanism for reporting the failure.

### 4. Termination is not rollback

Contract enforcement may terminate execution. Do not assume stack unwinding, RAII destruction or restoration occurs on contract termination. A violation handler may supply diagnostics under supported semantics, but it is not the transaction manager.

Supported recoverable failures use ordinary error returns or exception propagation according to ADR-0006 (ordinary error returns; only `std::bad_alloc` may propagate). RAII cleanup must correctly handle those paths. Throwing from a contract handler is not a general public error strategy, particularly across non-throwing functions or binary boundaries.

Mode-specific termination requirements are:

- **Ignored mode:** no contract termination occurs; all public validation, mandatory publication gates and supported rollback guarantees must still hold.
- **Observing mode:** observation must not be treated as rollback, recovery or guaranteed termination. The implementation must not publish a result that depends on the observer stopping execution.
- **Enforcing mode:** termination may occur for a contract violation. The process-termination path is outside the recoverable failure guarantee; no requirement may depend on destructors, rollback or post-termination cleanup running.

Optional notifications after publication cannot turn a committed result into an error that promises unchanged state. Process termination remains outside the existing strong failure guarantee; this does not relax any pre-publication safety requirement.

### 5. Constant evaluation verifies concrete cases

Make pure structural predicates, handle utilities and bounded transformations `constexpr` wherever practical within the selected toolchain and representation. Use `static_assert` for small independently specified fixtures: valid loop closure, edge splitting, restricted split/join equivalence, and rejection of wrong-loop or stale-handle inputs where representable at compile time.

`constexpr` and `static_assert` predate C++26; C++26 extends constant-evaluation capabilities and diagnostic facilities. Adoption does not require every transaction, allocation or numerical algorithm to become constant-evaluable.

Compile-time fixture success is evidence for those fixtures, not proof for arbitrary runtime models. Runtime failure injection, corruption tests, numerical adversarial tests and independent structural validation remain necessary. Large model fixtures should not be moved to compilation without measuring build-time and memory costs.

### 6. Reflection supports mechanical code generation

Reflection may generate member names, debug inspection, field enumeration and mechanical diagnostic formatting. Kernel semantics such as reciprocal incidence, radial membership, manifoldness, trim consistency and geometric partitioning remain explicitly specified predicates.

Do not derive stable public error identifiers, persistent schema field IDs or model serialization compatibility from incidental member order or source spelling. Such identities remain explicitly versioned. Reflection-generated checks must not silently omit a newly introduced semantic obligation merely because its storage representation is discoverable.

Reflection is a mandatory reference-toolchain language feature under ADR-0016. Using it to generate a particular validator or diagnostic is optional: an explicit implementation remains valid and must preserve the same observable behavior in every contract mode.

### 7. Toolchain and checking policy are explicit

Preserve ADR-0016: contracts, reflection and the other named language features are mandatory on the reference compiler; non-supporting compiler lanes remain non-blocking. Use its macros, feature probes and documented module/compiler workarounds rather than introducing an alternative contract abstraction. Feature-test macros and focused compile/run probes supplement the selected language-standard flag. This ADR neither changes the repository’s compiler baseline nor claims that all C++26 features are uniformly implemented.

Record the selected contract mode for each supported build profile: ignored, observing or enforcing. Apply a consistent intended policy across translation units and modules, while allowing deliberate diagnostic builds. Never infer contract behavior solely from `NDEBUG` or the name “release.”

Required safety and publication tests must pass with supplemental contracts ignored, observing and enforcing where the toolchain supports those modes. Tests must demonstrate that ignored mode still performs all mandatory validation and that observing mode does not provide the only failure path. Enabled-contract tests intentionally causing an enforcing violation run in isolated processes where supported. A portability fallback preserves ordinary checks and semantic guarantees, even when native annotation syntax is unavailable.

## Options

### Use native contracts for every public precondition

Rejected. Recoverable invalid input needs structured errors, and portable native-contract syntax alone does not guarantee non-ignorable enforcement or rollback.

### Use only ordinary checks and omit language contracts

Viable as a compatibility path, but not preferred when native contracts provide clear internal documentation and useful diagnostics without changing public behavior.

### Use native postconditions as the publication gate

Rejected when the body commits first. Validation must execute while edits are still private or reversible.

### Treat observing mode as sufficient enforcement

Rejected. Observation does not necessarily terminate execution, roll back state or prevent publication. Mandatory guarantees require explicit checks independent of observation.

### Treat enforcing mode as the public error mechanism

Rejected. Enforcement may terminate the process and is not a portable recoverable error policy. Public invalid input must be handled before any contract violation is required.

### Generate semantic validation entirely through reflection

Rejected. Discovering fields does not establish their geometric or topological meaning. Reflection-generated code may implement explicitly supplied rules, but cannot replace their specification.

### Require all kernel correctness checks at compile time

Rejected. User models and numerical inputs arrive at runtime. Constant evaluation complements runtime assurance and does not replace it.

## Consequences

- Contract notation improves visibility of internal obligations without making recovery depend on compiler settings or contract mode.
- Public validation, transaction and publication semantics remain stable across ignored, observing and enforcing configurations.
- Shared pure predicates reduce semantic drift between explicit checks, assertions and compile-time fixtures.
- Contract-mode testing must verify both diagnostic behavior and the continued independence of mandatory guarantees.
- Compiler feature gating and additional CI configurations require maintenance.
- Some checks may be repeated in diagnostic builds; optimize only after measurement while preserving mandatory gates.
- Reflection and constant evaluation are adopted where they reduce errors or maintenance cost, with no automatic runtime speedup assumed.

## Verification

1. Invalid public handles and unsupported requests produce the documented errors with contracts ignored, observing and enforcing where supported; none of these configurations mutates the input.
2. An injected mandatory-postcondition failure blocks publication even when native postconditions are ignored, not observed, or otherwise do not execute.
3. Recoverable failures during edits roll back or discard private state through supported return/exception paths, independently of contract handlers and termination behavior.
4. An intentional enforcing contract violation is diagnosed in an isolated-process test; the test does not assume destructor execution or rollback after termination.
5. Observing-mode tests confirm that a reported contract violation is not treated as a rollback signal or as the ordinary public error path.
6. Representative pure topology fixtures compile under `static_assert`, and corresponding runtime fixtures establish the same semantic results.
7. Predicate review confirms no correctness-critical side effects, invalid-handle dereferences or dependence on a particular number of evaluations.
8. Reflection-generated diagnostics preserve explicit stable identifiers; disabling reflection changes neither validation meaning nor public failure categories.
9. CI records compiler feature support and contract mode for each supported configuration. Merely compiling with a C++26 flag is insufficient evidence.

**Acceptance owner:** kernel API and toolchain maintainers, with topology-validation review.

## Integration

This proposal supplements accepted ADR-0016 and ADR-0006 without changing the mandatory language features, compiler baseline, contract macros or error model. Architectural acceptance remains human-owned under ADR-0015.

## References

- [C++ working draft: function contract specifiers](https://eel.is/c++draft/dcl.contract.func) — preconditions, postconditions and normal-return semantics.
- [C++ working draft: contract assertion evaluation](https://eel.is/c++draft/basic.contract.eval) — evaluation policies, violations, side effects and termination.
- [WG21 P2996R13: Reflection for C++26](https://www.open-std.org/jtc1/SC22/wg21/docs/papers/2025/p2996r13.html) — reflection mechanisms; diagnostic uses in this ADR are Stereon design choices.
- [GCC C++ standards support](https://gcc.gnu.org/projects/cxx-status.html) — feature-specific implementation status and feature-test information.

Sources were consulted in the preceding discussion on 2026-10-20. Working-draft and implementation-status URLs are mutable; this ADR intentionally avoids fixing compiler-version claims from them.

## Related decisions

- [ADR-0006](0006-error-handling.md)
- [ADR-0016](0016-cpp26-adoption-and-contracts.md)
- [ADR-0018](0018-validity-checking.md)
- [ADR-0039](0039-executable-topology-contracts.md)
- [ADR-0040](0040-atomic-topology-edit-publication.md)
- [ADR-0041](0041-layered-validation-evidence.md)
