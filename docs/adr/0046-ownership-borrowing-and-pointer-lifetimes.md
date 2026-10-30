# ADR-0046: Explicit ownership, borrowing and pointer lifetimes

- **Status:** Proposed
- **Date:** 2026-10-28
- **Owner:** Onur Tuncer
- **Phase:** 0, applying to all subsequent phases
- **Scope:** Kernel, adapters, bindings, tools and tests; third-party source remains governed by its own interface contracts
- **Related:** [ADR-0003](0003-geometry-polymorphism.md), [ADR-0004](0004-topology-storage-and-handles.md), [ADR-0005](0005-immutability-and-sharing.md), [ADR-0006](0006-error-handling.md), [ADR-0007](0007-threading-and-determinism.md), [ADR-0016](0016-cpp26-adoption-and-contracts.md), [ADR-0020](0020-cpp-c-abi-and-python-bindings.md), [ADR-0023](0023-logging-and-journal-sinks.md), [ADR-0038](0038-topology-traversal-lifetime-and-validation.md)

## Context

Stereon deliberately uses raw pointers for borrowed geometry, iterators, synchronous worker callbacks and type-erased storage. Ownership normally lives in values, containers, immutable snapshot tokens or smart pointers. A raw pointer is not inherently an error, and replacing every observer with shared ownership would add costs without addressing invalidation or data races.

The rules nevertheless need to be consistent across libraries. A pointer returned by `Target<T>()` can outlive an inline wrapper; a callback can leak an arena address despite a reference-return check; a diagnostic function can return a pointer into a shared string that another call replaces. `std::span`, references, iterators and `std::string_view` have the same lifetime issue. Const access does not retain storage or prevent its owner from changing it.

This ADR defines how ownership and borrowing are expressed, reviewed and verified. It complements ADR-0038's topology-specific lifetime rules; it does not weaken them or claim that existing code already satisfies every requirement.

## Options

1. **Ban raw pointers and use shared ownership everywhere.** Makes some lifetimes longer, but adds ownership traffic, permits cycles, and does not establish bounds, immutability or synchronization.
2. **Keep conventions local to each library.** Cheap initially; leaves temporary borrowing, callback retention and invalidation inconsistent and hard to review.
3. **Use explicit RAII ownership with scoped, documented borrowing.** Keep cheap observers where justified, retain ownership only across lifetime boundaries, and test invalidation and concurrency contracts.

## Decision

Adopt option 3. Apply these rules to new and changed interfaces; migrate existing violations through the inventory below. Acceptance requires the verification evidence, not merely publication of this document.

### 1. Express ownership separately from access

| Intent | Representation and rule |
| --- | --- |
| Owned data with local/value semantics | Value or owning container; prefer this before indirect allocation |
| Exclusive dynamic ownership | `std::unique_ptr<T>` or a resource-specific RAII type |
| Actual shared lifetime | `std::shared_ptr<T>`; prefer immutable payloads where the design requires immutability |
| Required borrowed object | `T&` or `const T&`; no ownership transfer |
| Optional borrowed object | `T*` or `const T*`; document what null means |
| Borrowed sequence/text | `std::span<T>` / `std::string_view`; document owner, bounds and invalidation |
| Topological identity | Typed handle plus the correct owning snapshot; a handle does not retain storage |
| Foreign/C ABI resource | Explicit opaque handle with documented create/retain/release rules, backed by RAII inside the adapter |

Raw C++ pointers do not transfer responsibility for freeing a dynamic allocation. Do not use a raw owning field when a standard or resource-specific RAII owner expresses the lifetime. A smart pointer parameter is appropriate when the function receives, transfers or retains ownership; an ordinary synchronous observer takes a reference, pointer or view instead.

Raw pointers inside RAII implementations, such as the model pointer into `SmallBuffer`'s own byte storage, are permitted. Their storage owner and object-lifetime protocol must be explicit. `const T*` and `shared_ptr<const T>` are not synchronization mechanisms; concurrent mutation through another alias still requires a valid synchronization or immutability contract.

### 2. Every escaping borrow has a lifetime contract

Document for every public borrowed return, stored observer field and callback interface:

- the object that owns the referenced storage;
- nullability and the valid element extent, where applicable;
- which operations invalidate the borrow, including destruction, assignment, move, mutation, reallocation and copy-on-write replacement;
- whether the callee may retain the borrow or pass it to another thread;
- what synchronization protects both access and lifetime.

Local pointers whose owner and use remain visibly within one scope need no repetitive comments. A stored observer must have a structural lifetime argument: its owner outlives it, or an accompanying ownership token retains the referenced storage. Ownership tokens must retain the actual pointee, not an unrelated object.

Do not rely on incidental survival, such as another shared owner happening to keep a payload alive, unless the interface promises that behavior. A vector's `reserve` is not a lifetime guarantee when later operations can exceed capacity. Copy-on-write sharing does not make references into mutable builder state stable.

### 3. Do not expose ephemeral storage through rvalue accessors

A member returning a borrow into its receiver must reject rvalue receivers unless the returned object independently retains ownership or the referenced storage is provably independent of the receiver's lifetime.

For example, geometry inspection is an lvalue operation:

```cpp
template <typename T>
[[nodiscard]] const T* Target() const & noexcept;

template <typename T>
const T* Target() const && = delete;
```

Apply the same principle to raw cache accessors, borrowed references, spans and views. It includes `GetModel()` if publicly exposed. Name an owning temporary before borrowing from it; where retaining the result is the real use case, provide an owning value or owning view instead. Deleting the rvalue overload also rejects immediate inspection of a temporary; callers migrate to a named owner.

Owning results remain available from temporaries. A retained topology range may outlive its originating Shape because it owns a snapshot token. Iterators still borrow their range under ADR-0038: no `enable_borrowed_range` opt-in, no iterator extraction from temporary ranges, and no iterator survival guarantee across range move, assignment or destruction. Ref-qualify range accessors where practical; range-based for over an owning temporary range remains supported.

Ref qualification prevents one class of dangling access, not all of it. Borrowing from a named object and then destroying or modifying that object remains caller misuse.

### 4. Edit callbacks cannot let mutable borrows escape

`PagedArena::Modify`, `ModifyPage` and similar controlled-edit callbacks grant access only for the callback invocation. Callbacks must not return, capture for later use, store elsewhere or dispatch asynchronously any pointer, reference, iterator or view into that storage. Returning an owning copy or unrelated value is permitted.

Keep re-entrant mutation/share guards active while the callback runs. Reject reference results and direct borrowed result forms, including raw pointers, spans, string views and reference wrappers, through a documented return-policy constraint. This conservative restriction may reject a pointer unrelated to the arena; such callers should return an owning value or an identifier instead.

A trait cannot discover every user-defined wrapper or captured address. Review and targeted analysis remain required; do not describe a return-type check as proving non-escape. Known wrappers carrying edit-storage borrows must not be approved as value results.

### 5. Callback and worker boundaries retain what they need

The low-level `void*` state and function-pointer bridge used by `ParallelFor` is permitted at the module/executor boundary. The implementation must complete or join all started work before returning or unwinding past borrowed state. The exact state type must match the callback cast. A null pointer is allowed only where its meaning is explicit, such as selecting serial execution with a null pool.

A stack work descriptor and pointers into caller-owned vectors are valid only while that completion guarantee holds. Allocate/size output storage before dispatch, prohibit structural vector changes during work, and give each task exclusive output elements or explicit synchronization. Separate logical elements do not necessarily mean separate memory locations: packed `vector<bool>` writes require special treatment.

Any future asynchronous interface must transfer or retain ownership of state and inputs until completion, including cancellation and error paths. Capturing a pointer or span by value does not retain its owner. Capturing a reference to a smart pointer does not retain an ownership copy either.

Mutable scratch is per invocation or per task/chunk, or explicitly caller-owned. Do not add shared mutable projector scratch or hidden caches to compensate for missing lifetime management. Respect ADR-0007's pool ownership, re-entrancy, nested-work and deterministic-order rules.

### 6. Encapsulate manual object lifetime and foreign interfaces

Placement construction, explicit destruction and pointer casts are permitted inside audited RAII/type-erasure or foreign-interface adapters. They must establish:

- sufficient size and alignment for the actual object;
- a pointer obtained through valid construction/access, with no use of stale pointers after reconstruction;
- type identity before recovering a concrete erased type;
- exactly-once destruction, valid move/self-assignment behavior and a safe state after construction failure;
- an owner that encloses the full lifetime of all internal observers.

Ordinary algorithm code must not introduce manual allocation/free pairs or cast away constness to mutate genuinely const data. Delegating a non-const accessor through its const overload is permitted when the underlying object is non-const. Character/byte views for I/O are permitted with correct extent and synchronous consumption or retained ownership.

Foreign functions may require raw owning handles, pointer/count pairs or matching allocation/free functions. Isolate that convention in an adapter; document allocator/deallocator pairing and ownership transfer on every success/failure path. Do not export C++ allocator ownership implicitly across a C ABI.

### 7. Diagnostics obey lifetime and concurrency rules too

A mutex that protects creation of a result does not protect an internal buffer after the function returns. General concurrent diagnostic APIs return an owning value or copy into caller-provided storage with a capacity and required-length contract. They must not return a mutable process-wide buffer as if it were stable.

The legacy debugger-only `stereon_describe_handle()` may temporarily retain its documented "until next call" behavior, explicitly restricted to serialized debugger invocation with no concurrent consumers. Add a caller-buffer alternative and migrate ordinary consumers to it. The legacy exception is named and scoped; it is not a model for new APIs or permission for thread-local result caches.

The snapshot registry may keep observers only while registration and deregistration coordinate with all lookups so storage cannot be destroyed during a lookup. Copy the result while the protection/retention is active; never let a borrow to a registry entry escape that protection accidentally.

### 8. Enforcement and errors

Borrowing preconditions for trusted C++ code are programmer contracts under ADR-0006 and ADR-0016. External null/length/handle data is checked at its boundary and reported through the applicable error API. No bounds or null check can prove that an arbitrary foreign address is live; foreign callers must meet the documented accessible-memory contract.

Use compile-time constraints/ref qualification where possible, checked lifetime/epoch instrumentation where justified, and native ASan/UBSan/TSan tests for runtime paths. Sanitizer success does not prove every lifetime safe. Do not add per-element shared ownership or full-model validation solely to avoid specifying a borrow contract.

New deviations require a named exception here or in a linked ADR with the owner, lifetime argument, affected interfaces and verification. Do not use local comments to silently weaken this decision.

## Migration inventory

The inventory is grounded in source revision `1776a52`. It is a work list, not evidence of completed enforcement. Its first disposition, with the evidence gathered and what stays open, is [the audit of 2026-10-28](../planning/adr-0046-audit-2026-10-28.md).

| Area | Required action |
| --- | --- |
| `Curve` / `Curve2d` / `Surface` inspection | Audit `Target<T>()`, `GetModel()` and related borrowed accessors; add lvalue-only borrowing or document an owning alternative; migrate temporary-inspection call sites |
| Trimmed/swept caches and Context pool accessor | Document owner/invalidation; prevent retaining a borrow from a temporary owner; keep the synchronous internal paths inexpensive |
| Topology ranges | Preserve ADR-0038 ownership tokens and non-borrowed-range semantics; enforce/test temporary-iterator restrictions |
| Arena edit callbacks | Extend return constraints to direct borrowed forms; document non-escape and review all callbacks |
| Parallel work descriptors | Verify join-on-completion/failure and pre-sized disjoint storage at every dispatch site |
| `SmallBuffer` and erased casts | Retain encapsulation; add/confirm alignment, lifetime, failure-path and move tests |
| Diagnostic registry and C string helper | Keep table lifetime protected; introduce caller-owned output and explicitly scope the legacy debugger exception |
| Bindings and foreign adapters | Document retained owner, copy policy and allocator/release pairing; no view may outlive its backing storage |

## Consequences

Borrowed access remains allocation-free and avoids reference-count work in hot loops. Ownership transfers become visible and escaping borrows become reviewable. This is not a raw-pointer ban or a mandate to put all storage on the stack.

Rvalue-access restrictions are source-breaking for some existing call sites; migrate them together and record the change under the pre-1.0 compatibility policy. Arena callback restrictions may also require callers to return values or typed IDs. The implementation and tests are follow-up work; this documentation change alone closes no memory-safety or release gate.

## Verification

| Criterion | Pass evidence |
| --- | --- |
| AC-0046-01: Interface inventory | Every public borrowed result, stored observer and callback boundary has an owner, invalidation and concurrency contract; exceptions are enumerated |
| AC-0046-02: Temporary borrowing | Compile-time tests accept lvalue inspection and owning temporary results, reject ephemeral rvalue inspection and temporary iterator extraction, and confirm topology ranges are not borrowed ranges |
| AC-0046-03: Retained lifetime | Native tests destroy/move original owners while retained ranges/owning results remain usable; documented invalidation is tested through compile-time or checked-instrumentation paths rather than silently dereferencing dangling objects |
| AC-0046-04: Edit non-escape | Compile-time tests reject reference, pointer, span and reference-wrapper callback results; tests cover callback failure and re-entrant edits; review covers capture/store escape paths that traits cannot prove |
| AC-0046-05: Parallel completion | Native tests verify no callback executes after return/unwind, owner release happens after completion, outputs remain disjoint, and 1/2/8/32-thread behavior meets ADR-0007; TSan reports no races in supported use |
| AC-0046-06: Erased storage | Native ASan/UBSan tests cover empty state, inline/shared payloads, supported alignment limits, copy/move/self-assignment, construction failure and exactly-once destruction |
| AC-0046-07: Concurrent diagnostics | Caller-buffer/owning results survive other diagnostic calls; concurrent lookup, snapshot release and output consumption pass lifetime tests and TSan; the legacy serialized-only helper is documented |
| AC-0046-08: Foreign boundaries | Tests cover null/length/handle validation, owner retention and matching release paths; contract tests do not treat arbitrary dangling foreign addresses as recoverable validated input |

**Acceptance:** keep this ADR Proposed until the owner reviews the policy, the migration inventory has disposition and the verification evidence is linked. New code can follow the proposed rules immediately; acceptance must not be inferred from the existence of this file.
