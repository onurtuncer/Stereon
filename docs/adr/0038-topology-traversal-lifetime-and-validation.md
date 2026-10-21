# ADR-0038: Topology traversal lifetime, invalidation and validation

- **Status:** Proposed
- **Revised:** 2026-10-20 (cross-reference executable contracts and validation evidence)
- **Date:** 2026-10-17
- **Owner:** Onur Tuncer
- **Phase:** Topology foundation
- **Related:** [ADR-0004](0004-topology-storage-and-handles.md), [ADR-0005](0005-immutability-and-sharing.md), [ADR-0006](0006-error-handling.md), [ADR-0007](0007-threading-and-determinism.md), [ADR-0016](0016-cpp26-adoption-and-contracts.md), [ADR-0018](0018-validity-checking.md), [ADR-0026](0026-modelling-operations-by-contracts.md), [ADR-0027](0027-controlled-topology-mutation-and-euler-operators.md), [ADR-0036](0036-typed-topology-ranges.md), [ADR-0037](0037-topology-traversal-semantics.md)

## Context

Typed IDs prevent mixing entity kinds but do not keep a store alive. An immutable range may outlive the Shape variable that created it; a builder cursor may survive rewiring even though every entity handle remains live. A malformed circular list can loop indefinitely. Entity generations alone do not solve these problems.

## Options

1. **Borrow every owner and rely entirely on caller discipline.** Cheap, but easy to dangle across temporaries, callbacks and bindings.
2. **Copy all topology into each traversal.** Safe but unacceptable for normal local queries.
3. **Retain immutable snapshots, keep iterators scoped to ranges and give private builder traversal a separate invalidation protocol.** Explicit and efficient without per-step ownership traffic.

## Decision

Adopt option 3.

### Frozen ranges

A public frozen range retains one immutable snapshot ownership token. The token refers to the complete immutable `TopologySnapshot` through its shared control block; it does not copy Shape’s page-pointer tables. The snapshot descriptor retains immutable pages, validation evidence and the logical topology store view.

Iterators borrow their range’s state. Destroying, moving or assigning the range invalidates its outstanding iterators; copying a range does not retarget those iterators. These types do not opt into `enable_borrowed_range`. Range-based for over a temporary range is supported through the lifetime of that range; extracting an iterator from a temporary is not supported. Document this distinction and test it.

Destroying or moving the original Shape variable does not invalidate a retained range. Editing a separate builder, `Freeze()` or `Compact()` producing another snapshot does not invalidate a range over the old snapshot. IDs extracted from a range still require the correct live owning snapshot for resolution. An iterator or ID is not a lifetime-owning substitute for its range/result.

Separate ranges and iterator copies can traverse one immutable snapshot concurrently. Sharing and mutating the same iterator object between threads is not supported. Maintained and prepared adjacency indexes are immutable after publication; iteration does not mutate shared caches or introduce global state.

### Builder cursors

Ordinary public topology ranges operate on frozen snapshots. Private mutation algorithms may use distinct builder cursor types only at structurally coherent points permitted by ADR-0027. They cannot expose partially edited transactions to ordinary consumers.

Each builder has a stable `BuilderToken` in addition to entity generations and a topology epoch. Any change to connectivity, membership, liveness, ordering or traversal-visible use orientation invalidates all its outstanding topology cursors before the change begins. Commit, rollback, Freeze, move and destruction also invalidate them. Rollback never restores an earlier epoch, even if it restores the topology. A failed operation that began a mutation invalidates cursors too.

Pure geometry or attribute changes may retain topology cursors only when their documented contract changes none of the fields those cursors observe. Conservative invalidation is permitted. Generation validity does not imply cursor validity.

Internal checked cursors verify owner token and epoch before accessing connectivity. A stale cursor is programmer misuse under ADR-0016, with a diagnostic naming the expected and actual owner/epoch. Do not present this check as making unsynchronised concurrent builder mutation safe. Builders remain single-threaded. Counter wrap must not allow stale epochs to match; choose a checked overflow strategy.

Do not support erase-while-iterating semantics initially. Collect IDs in a read phase, destroy the cursor, then mutate through the controlled layer; revalidate each planned target because an earlier edit may remove it. Such ID lists do not guarantee continuing adjacency or identity after edits.

### Validation and errors

Ordinary forward ranges assume structurally certified topology in the retained immutable snapshot. Every public query acquisition returns `std::expected<Range, KernelError>`. Raw typed-handle/store misuse inside the low-level C++ interface remains a contract violation as specified by ADR-0006; external numeric handles and imported data are checked at their boundary.

`Freeze()` is the authoritative validation gate. Import and deserialization must stage data and invoke the same validation before producing a frozen Shape. Builder commit establishes a coherent mutable state but does not establish ordinary traversal certification.

The snapshot stores a `ValidationCertificate` recording successful Freeze validation. Any topology mutation invalidates pending certification before the mutation begins. A new frozen snapshot receives a new certificate. Move transfers the certificate with the snapshot; serialization requires validation after deserialization rather than trusting a serialized certificate.

Do not rescan and validate the entire model on every local query. A valid immutable snapshot retains evidence of structural validation. Untrusted/deserialised topology must pass Freeze validation before it can supply ordinary ranges.

The validator uses a separate bounded diagnostic walk; it cannot assume the ordinary range is valid to prove its own validity. Before dereferencing links, check bounds, liveness, relation membership and required reciprocal links. Track revisits and use explicit finite bounds derived from store size to catch non-anchor cycles, missing closure and broken radial membership. Validation failure returns a structured `KernelError` with the offending relation/entity, never successful early termination.

Internal checked traversal may also detect corruption, but an ordinary range is not a fallible stream whose failure is encoded as end-of-range. If bounded inspection of incomplete topology is needed, provide `DiagnosticTopologyCursor` with separate `Item`, `End` and `Error` outcomes. It must not claim the ordinary forward-range contract or report a partial prefix as a complete query result.

## Resolved supporting decisions

### Snapshot ownership decision

A frozen Shape owns one immutable `TopologySnapshot` object. The snapshot contains the immutable logical store view, validation evidence, topology epoch identity and references to immutable page storage. Public ranges, recursive results, prepared indexes and binding cursors retain a shared ownership token for that snapshot.

The token is implemented as an intrusive or equivalent reference-counted control block selected by the topology module. It is not a copy of Shape’s page-pointer tables. Copying a range or result increments ownership only for the shared snapshot control block; iterator increment and dereference perform no ownership operations. A range retains the snapshot directly through the token, rather than borrowing the originating Shape variable.

The snapshot remains alive until the last Shape, range, recursive result, prepared index or binding cursor releases its token. Moving or destroying a Shape does not invalidate retained queries. A new Freeze or Compact creates a new snapshot and does not mutate or invalidate older snapshots. IDs and iterators do not own snapshots independently and cannot be resolved without the owning range, result or explicit snapshot handle.

The existing Arena representation is not assumed to provide this token. If it cannot provide shared immutable ownership without page-table copying, the topology module must add `TopologySnapshot` and its control block before public ranges are implemented.

Range construction must retain one token and must not perform O(number-of-pages) ownership work. The implementation must measure token construction, copying, retention and destruction costs. A range that retains a snapshot retains all storage reachable from that snapshot; callers are responsible for releasing long-lived ranges and results.

### Validation gate decision

`Freeze()` is the authoritative certification gate for ordinary traversal. A builder may not publish a frozen Shape until the topology validator has completed successfully. Imported or deserialized topology must enter a builder or equivalent staging object and pass the same validation before it can produce a frozen snapshot. Builder commit creates a coherent mutable state but does not certify it for ordinary public traversal.

The mandatory Freeze validation checks all structural invariants required by every ordinary topology query:

- bounds, liveness and generation validity of every referenced entity;
- ownership and relation membership of every link;
- reciprocal `Next`/`Prev` relationships;
- finite, anchored loop closure;
- finite, anchored radial-ring closure;
- required membership of coedges in loops and radial rings;
- required reciprocal links between topology records;
- root and disconnected-component membership;
- absence of references to tombstones or entities from another snapshot;
- validity of stored connectivity and traversal-visible orientation fields.

Category-specific invariants that are not required by every query remain deferred until the corresponding query is acquired. Such deferred checks must be named in that query’s contract and must return `KernelError`; they may not be silently treated as empty results.

Validation evidence is stored in the immutable snapshot descriptor as a `ValidationCertificate`. The certificate records the validator version, invariant set, source operation (`Freeze`), snapshot identity and successful completion. Any topology mutation invalidates the builder’s pending certificate. A new frozen snapshot receives a new certificate. Move operations transfer the certificate with the snapshot; they do not recreate or weaken it. Serialization does not transfer certification blindly: deserialization must validate the imported data before creating a certified snapshot.

The validator is independent of ordinary traversal. It uses bounded diagnostic walks with explicit visited state and finite limits derived from the relevant store sizes. Every relation type has a declared bound: loop and radial walks may visit at most the number of live coedges in the relevant store, and global validation may use the total live relation count plus one as the overflow sentinel. Counters use checked arithmetic and fail with a diagnostic error before overflow.

Ordinary forward ranges are acquired only from certified snapshots. They do not rescan the model. Query acquisition returns `std::expected<Range, KernelError>` for every public query, including queries whose inputs are typed handles, because snapshot certification, category-specific preconditions and preparation availability are part of the public contract. Low-level misuse of typed handles remains programmer misuse under ADR-0006; external numeric handles and imported data are validated at their boundaries.

Corruption is never represented as ordinary end-of-range. The separate diagnostic cursor reports `Item`, `End` or `Error`, and its error includes the relation, entity, snapshot identity and violated invariant. A partial diagnostic prefix is never reported as a complete successful query.

### Builder cursor ownership decision

Each builder owns a stable `BuilderToken` created at builder construction. The token is transferred by move and remains stable across builder moves; the moved-from builder is left without a valid token. Private builder cursors retain the token and the topology epoch observed at acquisition.

Every traversal-visible mutation invalidates cursors before the mutation begins, including failed mutations, rollback, commit, Freeze, move and destruction. Rollback never restores an earlier epoch. Entity generations validate entity handles only; they do not validate builder cursors.

Checked cursor access compares both the stable builder token and the captured topology epoch. A mismatch is programmer misuse under ADR-0016 and produces a diagnostic naming the expected and actual owner/epoch. Builders remain single-threaded; this check does not make concurrent mutation safe.

## Acceptance criteria

### AC-0038-01 — Snapshot retention after Shape destruction

- **Requirement:** A retained range, recursive result, prepared index or binding cursor keeps its immutable snapshot alive after the originating Shape is moved or destroyed.
- **Complexity:** Retention and release perform O(1) control-block ownership work with respect to page count.
- **Owner:** Snapshot ownership implementation owner.
- **Validation:** Create a Shape, acquire each retained-query type, destroy or move the Shape, then traverse or resolve through the retained object under sanitizers. Verify exact sequence and snapshot identity. Instrument page-table copies and ownership operations; reject any O(number-of-pages) retention path.

### AC-0038-02 — Snapshot isolation across Freeze and Compact

- **Requirement:** Freeze or Compact creates a new snapshot without mutating or invalidating older retained queries.
- **Owner:** Snapshot and builder integration owner.
- **Validation:** Acquire a range from snapshot A, mutate a builder, Freeze or Compact to snapshot B, then traverse A and B. Verify A’s sequence, certificate and identity remain unchanged and B reflects only its own state. Compare entity IDs across compaction only through the explicit remap.

### AC-0038-03 — Range and iterator lifetime

- **Requirement:** Iterators borrow their range state; range destruction, move and assignment invalidate outstanding iterators, while range copies retain independent traversal positions.
- **Owner:** Range implementation owner.
- **Validation:** Add compile-time checks that ranges do not model `borrowed_range`. Runtime-test temporary range use in range-based for, iterator extraction from temporaries, range move/assignment, copied ranges and copied iterators. Invalid use must be rejected by the documented contract or detected under debug instrumentation; no use-after-free is permitted under sanitizers.

### AC-0038-04 — Builder cursor invalidation

- **Requirement:** Every traversal-visible mutation invalidates builder cursors before mutation begins, including failed mutations, rollback, commit, Freeze, move and destruction.
- **Complexity:** Checked cursor access performs O(1) owner-token and epoch comparisons before following connectivity.
- **Owner:** Builder and cursor implementation owner.
- **Validation:** Save cursors, then exercise insert, delete, relink, reanchor, reverse, commit, rollback, Freeze, move, destruction and failed mutation paths. Verify stale access is detected before any stale link is read. Specifically verify rollback to byte-for-byte equivalent topology still rejects the old cursor.

### AC-0038-05 — Token and epoch identity

- **Requirement:** Builder moves transfer the stable `BuilderToken`; moved-from builders have no valid token; epoch counters cannot wrap into a false match.
- **Owner:** Builder ownership implementation owner.
- **Validation:** Test move construction, move assignment, moved-from access, repeated mutations and counter-boundary behavior using a checked test configuration. Verify diagnostics identify expected and actual owner/epoch values.

### AC-0038-06 — Freeze certification gate

- **Requirement:** No ordinary public range or result can be acquired from an uncertified, imported or deserialized topology.
- **Owner:** Validation and import integration owner.
- **Validation:** Attempt queries before Freeze, after builder commit but before Freeze, after failed Freeze, from malformed imports and from deserialized data carrying a forged certificate. Verify `std::expected` contains `KernelError`, no ordinary range is returned and the error identifies the failed certification or precondition.

### AC-0038-07 — Bounded independent validation

- **Requirement:** `ValidateTopology(snapshot)` and diagnostic cursors terminate within declared finite bounds for every relation type and never rely on ordinary traversal.
- **Complexity:** Loop and radial walks visit at most the relevant live-coedge count; global validation uses the declared live-relation bound plus one as overflow sentinel. Counters use checked arithmetic.
- **Owner:** Validator implementation owner.
- **Validation:** Corrupt `Next`, `Prev`, radial links, owner IDs, tombstone references, generation values, reciprocal links and root membership independently. Include non-anchor cycles, missing closure and oversized counters. Verify termination within the bound, no hang, no out-of-bounds access and structured `KernelError` naming relation, entity, snapshot identity and invariant.

### AC-0038-08 — Diagnostic cursor state machine

- **Requirement:** `DiagnosticTopologyCursor::Next()` distinguishes `Item`, `End` and `Error`; corruption is never reported as `End`, and a partial prefix is never reported as complete success.
- **Owner:** Diagnostic API owner.
- **Validation:** Inject corruption before the first item, between valid items and at closure. Verify exact state transitions and error payloads. Confirm ordinary forward ranges never expose this fallible stream behavior.

### AC-0038-09 — Concurrent immutable traversal

- **Requirement:** Independent ranges, iterator copies and prepared indexes may be traversed concurrently after publication without races or sequence variation. Sharing one mutable iterator object remains unsupported.
- **Owner:** Threading and topology query owners.
- **Validation:** Run repeated concurrent traversals under thread sanitizers with multiple ranges, copied ranges and prepared indexes. Verify identical sequences, no races, no mutation of shared caches and no hidden global state.

### AC-0038-10 — Retained-memory accounting

- **Requirement:** Snapshot retention and optional prepared-index retention are measurable and reclaimable according to the ownership policy.
- **Owner:** Snapshot and performance owners.
- **Validation:** Measure retained bytes for ranges, results, preparations and binding cursors across page counts. Verify local range construction does not copy every page-pointer table. Release owners in permutations and confirm snapshot and optional indexes are reclaimed only after the final required owner releases them.

### AC-0038-11 — Integration review gate

- **Requirement:** Implementation measurements, ownership traces, invalidation tests, validation diagnostics and import tests are attached to the ADR review.
- **Owner:** Topology integration owner.
- **Validation:** Keep ADR-0038 Proposed until all acceptance criteria pass and reviewers confirm compliance with snapshot, builder, validation and diagnostic API decisions.

## Consequences

Frozen queries remain valid while later edits produce new snapshots. Retaining old ranges retains old snapshot memory, which callers must release deliberately. Range moves have explicit iterator consequences. A lightweight snapshot descriptor and builder epoch infrastructure are required, and validation walkers must remain independent of trusted traversal.

This policy extends existing ownership and contract decisions without changing the four-byte release ID representation. It does not make IDs globally unique or introduce stable iterators across builder edits.


## Related executable contract proposals

- [ADR-0039: Make topology-operator contracts executable and auditable](0039-executable-topology-contracts.md)
- [ADR-0040: Stage topology edits and publish complete results atomically](0040-atomic-topology-edit-publication.md)
- [ADR-0041: Use layered validation with revision-bound evidence](0041-layered-validation-evidence.md)
- [ADR-0042: Assign C++26 contracts and kernel validation distinct responsibilities](0042-cpp26-contract-validation-responsibilities.md)
