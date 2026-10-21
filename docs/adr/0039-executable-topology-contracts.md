# ADR-0039: Make topology-operator contracts executable and auditable

- **Status:** Proposed
- **Date:** 2026-10-20
- **Owner:** Onur Tuncer
- **Phase:** 0
- **Related:** Operation contracts; topology/geometry separation; controlled mutation; accepted operator catalogue in ADR-0027 §9.

## Context

A mathematical Euler count does not establish correct incidence, ownership or geometric embedding. A function can add exactly one edge and one face while leaving a broken loop or updating the wrong face. Prose requirements also become ambiguous when callers cannot tell whether they are checked, assumed or established by previous validation.

The existing operation-contract proposal requires documented preconditions, postconditions and failure behavior. This ADR specifies how those obligations are implemented without requiring a full geometric revalidation before every local edit.

## Decision

Every supported topology operator shall have a stable contract identifier, a supported-domain profile and a reviewed enforcement specification. Every normative requirement shall name its enforcement mechanism: direct check, valid prior evidence, or preservation by a reviewed construction rule. An unchecked caller obligation must be explicitly identified and may not be advertised as verified.

### Contract contents

Each operator specifies:

1. Inputs, model/builder identity, accepted topology category and supported cases.
2. The read neighborhood and expected write footprint, including dependent associations, incidence indexes and metadata.
3. Preconditions with stable identifiers and failure classifications.
4. Success postconditions, including correspondence to the input and required lineage.
5. Frame conditions: which observable properties and entities remain unchanged.
6. Invariants temporarily suspended inside its private edit scope.
7. Evidence required or produced, including its scope and invalidation rules.
8. Failure behavior, builder reusability, and handle/cursor invalidation.

This specification may live in code-adjacent documentation and tests; a generated contract DSL is not required.

### Public checks and trusted internal paths

Public Result-returning operations and queries enforce handle ownership, liveness/generation, parameter validity and supported-profile selection before unsafe access or mutation. Low-level value-returning accessors and Context modifiers retain the programmer-obligation contract tier of ADR-0006 and ADR-0016; this proposal does not add a Result or an always-on branch to every field read. Assertions disabled in release builds shall not be the sole enforcement mechanism for these checks.

Expensive conditions may be established by evidence valid for the exact input revision and numerical context. A trusted internal path is private or explicitly designated; it accepts validated inputs or a prepared plan bound to those inputs. It preserves the handle-integrity and checked-build requirements of ADR-0004, ADR-0006 and ADR-0016. A boolean such as `assume_valid=true` without a defined evidence/lifetime policy is insufficient.

No unrestricted setters for authoritative Next/Prev/Radial, ownership or maintained incidence are exposed to modelling algorithms. Controlled mutation remains the governing boundary.

### Preconditions, postconditions and errors

A violated input requirement yields the existing public error mechanism with an appropriate semantic category: invalid input, unsupported profile, unresolved numerical condition, cancellation/resource failure, or an invalid result detected by explicit validation. This ADR does not introduce a competing error hierarchy.

An explicit result-validation failure returns `KernelError::InvalidResult` with its report under ADR-0006 and ADR-0018 and is counted as a kernel bug. A numerical verification stage may instead report unresolved according to its contract. Native contract violations remain programming errors under ADR-0016 and are never converted into `KernelError`. No failed mandatory result check may publish a nominally successful shape.

Diagnostics identify the operation, contract condition, stage and relevant entities when available. Diagnostic detail must not require dereferencing an invalid handle. Diagnostics for discarded candidate entities identify them as candidate-local rather than reusable output handles.

### Example: restricted SplitDiskFace

This example defines a restricted MEF-style geometric wrapper; it does not replace or narrow the accepted classical Euler catalogue in ADR-0027 §9. The illustrative name SplitDiskFace is not a claim of an implemented API. The wrapper accepts one regular disk face, one simple exterior loop with no holes, two distinct nonadjacent boundary occurrences, and a supplied simple cut. The cut is certified to meet the existing boundary only at its endpoints, with its interior inside the face. Seam-crossing and singular configurations require separate supported profiles. The topology core consumes the prepared cut specification; it does not solve the intersection problem.

| ID | Requirement | Establishment |
| --- | --- | --- |
| PRE-01 | All handles are live and belong to the selected builder/revision. | Checked before dereference. |
| PRE-02 | The selected occurrences lie on the same simple exterior loop and are distinct and nonadjacent. | Structural check or matching scoped evidence. |
| PRE-03 | The face and cut satisfy the declared disk and geometric profile. | Geometric preparation evidence under an explicit numerical context. |
| PRE-04 | Attribute, lineage and association policies permit the split. | Policy evaluation in the enclosing transaction. |
| POST-01 | Net changes are ΔV=0, ΔE=1, ΔFace=1, ΔLoop=1, ΔCoedge=2. | Count comparison for this profile only. |
| POST-02 | Two closed directed boundary cycles replace the old one; every old boundary occurrence belongs to the correct child exactly once. | Local connectivity/correspondence checks. |
| POST-03 | The divider has exactly the two intended opposite coedge uses; unaffected radial relationships remain correct. | Reciprocal incidence checks. |
| POST-04 | Children inherit the required surface, orientation and shell membership; associations and trims agree. | Association checks and prepared construction evidence. |
| POST-05 | Child regions partition the original region, with disjoint interiors, within the declared assurance. | Reviewed geometric construction and/or independent semantic validation. |
| POST-06 | Lineage and attributes are complete; frame conditions hold. | Transaction result and policy checks. |

Net entity counts do not decide whether a parent identity is retained or retired; the existing identity/lineage policy controls that choice. Boundary occurrences, not vertex identifiers alone, determine insertion positions. Required lineage may describe one parent yielding two children without promising identical storage slots.

### Performance and constant evaluation

Mandatory data checks and publication gates remain active in release builds. Per-operator local topology checks run in debug and sanitizer builds as required by ADR-0027 §9; release relies on the independent Freeze/Commit validator. Contract-specific pre-publication obligations remain mandatory and may be established by reviewed construction or current evidence where permitted. Validation cost is measured against affected coedges and dependencies, not assumed constant merely because an operator is called local. Pure predicates and bounded topology transformations should be `constexpr` where their dependencies and supported language/toolchain permit it. Runtime allocation, evidence storage and transactions need not become constant-evaluable to satisfy this ADR. Constant evaluation does not weaken contracts.

## Options

- **Documentation only:** rejected for directly checkable public safety requirements.
- **Debug assertions only:** rejected because production input errors and structural defects would escape the advertised boundary.
- **Full-model validation at every primitive:** rejected as the default; publication-level validation remains governed by ADR-0041.
- **Require a language contracts feature:** deferred; the semantics must work independently of that facility.

## Consequences

Contracts become traceable to code and tests, and internal optimizations can reuse valid evidence. The costs are contract maintenance, structured diagnostics and explicit dependency bookkeeping. A reviewed construction argument is useful evidence but is not a claim of machine-checked formal proof.

## Verification

- Every supported operator has the required contract fields and an enforcement entry for each requirement.
- Invalid handles, wrong-builder handles and unsupported profiles fail before mutation or unsafe access in release builds.
- SplitDiskFace fixtures independently specify boundary cycles, opposite uses and frame conditions; count-only verification is insufficient.
- Tests detect a wrong loop splice even when entity deltas are correct.
- A concave-face cut leaving the face is rejected or reported unresolved by the geometric wrapper despite plausible topology.
- Split followed by its restricted inverse restores semantic structure under an explicit entity correspondence; tests do not require reused storage indices.

**Acceptance owner:** topology operator maintainer, with independent validator review.

## Related decisions

- [ADR-0006](0006-error-handling.md)
- [ADR-0016](0016-cpp26-adoption-and-contracts.md)
- [ADR-0018](0018-validity-checking.md)
- [ADR-0026](0026-modelling-operations-by-contracts.md)
- [ADR-0027](0027-controlled-topology-mutation-and-euler-operators.md)
- [ADR-0038](0038-topology-traversal-lifetime-and-validation.md)
- [ADR-0040](0040-atomic-topology-edit-publication.md)
- [ADR-0041](0041-layered-validation-evidence.md)
- [ADR-0042](0042-cpp26-contract-validation-responsibilities.md)
