# ADR-0040: Stage topology edits and publish complete results atomically

- **Status:** Proposed
- **Date:** 2026-10-20
- **Owner:** Onur Tuncer
- **Phase:** 0
- **Related:** Strong failure guarantee in operation contracts; immutable snapshots; controlled mutation; iterator epochs; lineage and attribute policies.

## Context

A topology edit touches more than edge links. Associations, roots, indexes, tolerances, history and attributes must remain mutually consistent. An allocation failure or failed validation after a successful splice must not leak a partial model. Reversing an Euler operator is not a sufficient general rollback mechanism: the inverse may require information already lost or preconditions that no longer hold.

## Decision

Supported modelling operations use an explicit private edit scope and one publication boundary. Use the candidate-copy transaction selected in accepted ADR-0027 §9: edit a private builder copy whose pages are shared until written, replace the original builder on success, and discard the candidate on failure. This ADR refines that lifecycle; it does not reopen storage selection or introduce an undo journal.

### Lifecycle

1. Validate input identity and establish the supported domain.
2. Prepare a plan and obtain required evidence bound to the input state.
3. Open a private edit scope; reserve resources where practical.
4. Apply coordinated changes and record affected entities/dependencies.
5. Establish required postconditions, attributes, lineage and validation evidence.
6. Freeze and publish the completed result through one controlled boundary.

For immutable inputs, plan evidence remains bound to that snapshot. For mutable builders, exclusive edit ownership or a revision check prevents changes between planning and execution. A stale plan is rejected or recomputed; it is not used opportunistically.

The actual publication step must not perform work that can expose half-published state. Required allocations and validation happen before publication; implementation-specific synchronization establishes a single commit point. Concurrent mutable writers are not introduced by this ADR. If unsupported, they are rejected by the builder ownership policy.

### State covered by atomicity

Atomicity covers topology, geometric associations, changed geometry, roots, maintained incidence, public entity identity, tolerance metadata, lineage and attribute results. Shared geometry visible through input snapshots is immutable or protected by the existing copy-on-write policy. An unchanged input must remain unchanged through every supported alias.

Diagnostics, allocation statistics and semantically transparent caches may change. A cache must remain correct or be invalidated; cache state must not change modelling semantics.

### Builders and cursors

Immutable published snapshots and their supported readers remain valid after a failed operation. A caller-reusable builder retains its entry logical model state on recoverable failure: the private candidate is discarded under ADR-0027 §9. This proposal does not add a failure mode that poisons the caller's previously valid builder. Public modelling operations still preserve their input models under the strong failure guarantee.

Before the first mutation, affected builder cursor epochs advance according to the iterator policy. Rollback does not decrease epochs or revive old cursors. This is an explicit builder-lifetime effect, not permission to mutate an immutable input snapshot. Handles created inside an aborted edit never become valid public handles; slot reuse must not resurrect them through generation collisions.

Failure recovery must not depend on fresh allocations or on a general-purpose inverse modelling operation. Discard the private candidate, preserving the original builder state under ADR-0027. Candidate cleanup must not throw. A failed private candidate is discarded and cannot be reused or published as a successful model.

Allocation failures follow ADR-0006: the kernel does not throw, but `std::bad_alloc` may propagate. Atomicity applies whether the failure is reported as an error result or by supported exception propagation. Process termination and hardware failure remain outside the guarantee.

### Composition, cancellation and callbacks

Inner Euler edits do not independently publish when participating in one outer modelling operation. Inner commits mean acceptance into the private candidate, not publication. An outer failure discards or restores all constituent edits. A separately committed multi-step user workflow must say so explicitly.

Cancellation observed before the commit point produces no published change. Once the commit point has passed, return the committed result; a late cancellation must not produce a failure response that falsely implies nothing changed.

No user callback or observer may see or retain authoritative handles to partial candidate state. Progress callbacks carry safe progress/diagnostic data only and may not reenter the same mutable builder. Optional notifications after publication cannot retroactively turn a committed operation into an unchanged-state failure.

### Illustrative control flow

```text
plan = validate_and_prepare(input, request)
edit = begin_private_edit(input, plan.revision)
candidate = apply_coordinated_edit(edit, plan)
establish_postconditions(candidate, plan)
resolve_lineage_attributes_and_required_evidence(candidate)
result = freeze_and_prepare_publication(candidate)
publish(result)  // the sole commit point
return success(result)
```

Every supported failure before publication leaves the public input unchanged. A concrete RAII implementation must guarantee rollback/discard during stack unwinding; this pseudocode does not supply that mechanism.

## Options

- **Mutate first and let callers repair:** rejected for public operations.
- **Use the inverse Euler operation as undo:** rejected as a universal strategy because metadata and identity can be lost.
- **Copy the entire model every time:** permitted initially for correctness but not mandated; benchmark the accepted copy-on-write candidate mechanism.
- **Publish each primitive separately:** allowed only in explicitly separate public operations, not inside an atomic composite.

## Consequences

Failures are recoverable and model snapshots remain reliable. Transactions require memory and careful cleanup engineering. Cursor invalidation survives rollback, so callers must reacquire builder cursors after a mutation attempt. Publication, cancellation and error reporting must agree on whether commit occurred.

## Verification

- Inject failures at allocation, splice, association update, attribute resolution, lineage creation, validation and pre-publication stages; public inputs remain semantically unchanged through aliases.
- Caller builders retain logical state; discarded private candidates cannot escape; stale cursors remain invalid after rollback.
- Aborted candidate handles remain invalid after later successful allocations.
- A composite edit that fails after several successful inner operations publishes none of them.
- Rollback succeeds under allocation failure without allocating recovery storage.
- Cancellation immediately before and after the commit point yields outcomes consistent with actual publication.

**Acceptance owner:** mutation/storage maintainer with iterator and public API review.

## Related decisions

- [ADR-0006](0006-error-handling.md)
- [ADR-0016](0016-cpp26-adoption-and-contracts.md)
- [ADR-0018](0018-validity-checking.md)
- [ADR-0026](0026-modelling-operations-by-contracts.md)
- [ADR-0027](0027-controlled-topology-mutation-and-euler-operators.md)
- [ADR-0038](0038-topology-traversal-lifetime-and-validation.md)
- [ADR-0039](0039-executable-topology-contracts.md)
- [ADR-0041](0041-layered-validation-evidence.md)
- [ADR-0042](0042-cpp26-contract-validation-responsibilities.md)
