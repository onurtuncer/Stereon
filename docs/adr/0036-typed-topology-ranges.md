# ADR-0036: Typed read-only topology ranges

- **Status:** Proposed
- **Date:** 2026-10-17
- **Owner:** Onur Tuncer
- **Phase:** Topology foundation; implement with the topology store before dependent modelling algorithms
- **Related:** [ADR-0004](0004-topology-storage-and-handles.md), [ADR-0005](0005-immutability-and-sharing.md), [ADR-0006](0006-error-handling.md), [ADR-0007](0007-threading-and-determinism.md), [ADR-0017](0017-naming-conventions.md), [ADR-0020](0020-cpp-c-abi-and-python-bindings.md), [ADR-0025](0025-separate-topology-from-geometry.md), [ADR-0037](0037-topology-traversal-semantics.md), [ADR-0038](0038-topology-traversal-lifetime-and-validation.md)

## Context

ADR-0004 calls for adjacency views such as `FacesOf(edge)`, `EdgesOf(vertex)` and loop walks. Algorithms need a common interface without depending on page layout, tombstones or raw connectivity fields. A single generic explorer would obscure whether it returns immediate children, all descendants, entity identities or oriented uses. Exposing vectors for every query would impose allocation even on simple local walks.

## Options

1. **Expose storage arrays and raw links.** Minimal wrappers, but clients duplicate traversal rules and depend on physical representation.
2. **Return materialised vectors for every query.** Easy ownership and random access, at the cost of eager work and allocation.
3. **Provide typed ranges for validated local traversal and explicit materialised results for expensive queries.** Keeps common walks light while exposing cost and semantics.

## Decision

Adopt option 3. Export traversal through the topology module, using the existing module/partition conventions. Do not place topology traversal in geometry, presentation or the low-level handle module.

The topology module owns maintained structural indexes. Explicit prepared immutable indexes are owned by `TopologyQueryPreparation`, which retains the snapshot token and is created only by an explicit preparation call. Ordinary local range construction never allocates a preparation object or performs a hidden whole-store scan.

### Query families

Use separate names for these contracts:

|Query                                        |Result and scope                                                                         |
|---------------------------------------------|-----------------------------------------------------------------------------------------|
|`AllEntities<K>(shape)`                      |Every live entity of kind K in this snapshot’s store, once; not a root-reachability query|
|`LoopsOf(shape, face)`                       |Immediate loop membership of the face                                                    |
|`CoedgesOf(shape, loop)`                     |One finite walk around the loop in its stored direction                                  |
|`RadialCoedgesOf(shape, edge)`               |All coedge uses of that edge, including repeated face uses                               |
|`FacesOf(shape, edge)`                       |Distinct incident faces, once each; loses use multiplicity deliberately                  |
|`EdgesOf(shape, vertex)`                     |Distinct incident edges; a closed edge is returned once                                  |
|`CollectDescendants<K>(shape, root, options)`|Explicit eager recursive query with declared root inclusion, deduplication and ordering  |

A collection over the whole store is not presented as descendants of a particular solid. A Shape may contain multiple roots and disconnected components; `AllEntities<K>` enumerates them all without implying containment.

Range elements are typed IDs or small immutable use values. Dereferencing never returns a mutable topology record. Cross-kind conversions remain ill-formed. Handle values are meaningful only with their owning snapshot; a collected vector of IDs is not a standalone model.

### Iterator contracts

1. Local immutable range types shall model `std::ranges::forward_range`: copied iterators have independent positions and repeated traversal yields the same sequence. Dereference may return a handle by value. Do not promise writable references, contiguous storage, random access or bidirectional traversal merely because some underlying fields support them.
2. Range construction, `begin`, increment and dereference must meet the complexity requirements of the concepts actually advertised. If finding the first live slot or preparing unique adjacency is nonconstant work, do it explicitly at construction or preparation and document it; do not hide a whole-model scan in `begin()`.
3. Implement `size()` only when exact size is available with the advertised complexity. Scanning a ring to count it is an explicit linear query, not a constant-time size member. Empty queries are valid ranges with `begin() == end()`.
4. Standard customization names retain their standard spelling (`begin`, `end`, `size`, iterator associated types); domain functions follow ADR-0017.
5. Do not implement fundamental local traversal through coroutine generators. Independent forward iteration, explicit ownership and predictable allocation are the initial requirements. A future single-pass streaming API may use coroutines under a distinct contract supported by measurements.
6. No implicit parallel traversal. Consumers may partition independent work and merge in the fixed order required by ADR-0007.

### Performance and storage boundary

For a table with S allocated slots and L live entities, an alive-bitset scan costs O(S) in the conservative bound, not O(L). Local loop and radial walks with k uses cost O(k) and require O(1) traversal state after snapshot retention. They must allocate no memory and perform no ownership-reference increments per step.

Distinct adjacency can need scratch space. Implement it through a prepared immutable index or an explicit eager result; disclose preparation time and memory. Do not claim an allocation-free `FacesOf` range if it builds a hidden hash set. Vertex incidence requires a maintained index or explicit preparation: repeatedly scanning every edge is not an acceptable hidden implementation of a local adjacency query.

Recursive discovery is an explicit operation returning `std::expected<OwnedTraversalResult<K>, KernelError>` in the sketch API. The result retains its snapshot and stores ordered IDs. For a unique reachable-subgraph traversal with V entities and A incidences, target O(V+A) time and O(V) visited state, excluding any documented sorting. Recursive occurrence enumeration has output-sensitive cost and may be much larger than unique traversal.

Internal page spans remain an optimization boundary, not public topology iteration. No public range exposes page size or mutable storage. Benchmark page-oriented implementations without changing the observable traversal contract.

## Resolved supporting decisions

### Topology ownership decision

Each frozen Shape owns one complete immutable logical topology snapshot. Topology stores, roots, disconnected components, entity tables, incidence relationships and cross-store relationships belong to that snapshot’s topology module. Topology stores are not shared mutably between Shapes. Structural sharing of immutable pages is permitted as an implementation optimization, but the logical store view and ownership boundary remain Shape-specific.

A Shape may contain zero, one or many independent topology roots and disconnected components. Roots are represented by explicit root records in the snapshot, grouped by topology category where required. Root enumeration is an explicit whole-snapshot query and is not inferred from descendant traversal. `AllEntities<K>(shape)` enumerates every live entity of kind K in the snapshot exactly once, regardless of root membership. It is therefore never described as a descendant query.

The topology module owns maintained vertex-incidence indexes and other structural adjacency indexes required by local queries. Prepared distinct-adjacency indexes are owned by an explicit `TopologyQueryPreparation` object associated with one snapshot. They are not stored in mutable Shape state and are not algorithm-specific global caches.

`TopologyQueryPreparation` retains the snapshot token, is immutable after construction and may be shared across threads after publication. It owns only indexes explicitly requested by the caller. Its lifetime ends when the last preparation object or result that embeds its prepared data is released. The snapshot may reclaim no prepared index while a preparation object or result still references it. Ordinary local ranges do not create or retain preparation objects implicitly.

Maintained indexes are built or updated by the topology builder and published with the frozen snapshot. Explicit prepared indexes are built from a frozen snapshot by the topology query-preparation API. Hidden whole-store scans, hidden hash sets and hidden allocation in ordinary local range construction are rejected because they obscure complexity, create unpredictable ownership traffic and make local queries depend on unrelated model size.

### Prepared-index ownership decision

Maintained structural indexes belong to the immutable topology snapshot and are published with it. Explicit distinct-adjacency preparation belongs to `TopologyQueryPreparation`, which retains the snapshot token and owns its immutable index data. Recursive traversal owns its visited state and ordered output in `OwnedTraversalResult<K>`; it does not install a cache in the snapshot.

Prepared indexes are never hidden inside ordinary ranges. They may be reclaimed when the last preparation object or result using them is destroyed. A snapshot does not retain optional prepared indexes merely because one query used them. This policy is chosen to prevent unbounded snapshot retention, avoid hidden global caches and make preparation cost and memory visible to callers.

### API naming decisions

The normalized names are accepted:

- `AllEntities<K>(shape)` for whole-snapshot entity enumeration;
- `LoopsOf(shape, face)`;
- `CoedgesOf(shape, loop)`;
- `RadialCoedgesOf(shape, edge)`;
- `FacesOf(shape, edge)`;
- `EdgesOf(shape, vertex)`;
- `CollectDescendants<K>(shape, root, options)`.

The `Of` convention is retained because it identifies the source relation without implying containment, geometric ordering or manifoldness. Entity queries return unique entity IDs. Use queries preserve incidence multiplicity.

The manifold-only helper is named `ManifoldOppositeFace`. It validates that the edge has exactly two valid coedge uses, that the uses belong to distinct faces and that the required manifold category invariants hold. It returns `std::expected<FaceId, KernelError>` and never exists as an unqualified `Twin` or `OppositeFace` operation.

The orientation helpers are named:

- `StoredEdgeUseSense(coedge)`;
- `TraversalReversal(use)`;
- `EffectiveEdgeDirection(use)`.

`EffectiveEdgeDirection` composes stored sense and traversal reversal exactly once. Callers must not manipulate packed orientation bits directly.

The diagnostic API is named `ValidateTopology(snapshot)` for complete bounded validation and `DiagnosticTopologyCursor` for bounded inspection. `DiagnosticTopologyCursor::Next()` returns a discriminated result with `Item`, `End` and `Error` states. These names are separate from ordinary range APIs.

C and Python adapters may choose binding-specific spellings, but they must preserve the native distinctions between entity IDs, oriented uses, occurrence context, ranges, preparation objects and snapshot ownership. Native C++ names are not required to appear unchanged in bindings.

## Acceptance criteria

The implementation shall not advance from Proposed until all criteria below are satisfied and recorded in the ADR implementation review.

### AC-0036-01 — Ownership and allocation

- **Requirement:** Every public local range retains exactly one snapshot ownership token through a shared control block and does not copy Shape page-pointer tables.
- **Complexity:** Range construction performs O(1) ownership operations with respect to the number of pages. Iterator increment and dereference perform O(1) traversal work for local loop and radial ranges and perform zero ownership increments/decrements.
- **Allocation:** Local loop and radial range construction, `begin`, increment and dereference perform zero dynamic allocations. Distinct adjacency must use either an explicitly prepared index or an eager result whose allocation is visible in the API.
- **Owner:** Topology range implementation owner.
- **Validation:** Instrument token operations, allocations and page-table copies. Test snapshots with 1, 10, 1,000 and 10,000 pages. Pass only if ownership-operation count is constant across page counts and no page-table copy occurs.

### AC-0036-02 — Range concepts and iterator lifetime

- **Requirement:** Each local range models exactly the advertised standard concepts and is non-writable.
- **Complexity:** `begin`, increment and dereference satisfy the documented bounds; no hidden whole-store scan occurs in `begin`.
- **Owner:** Topology range implementation owner.
- **Validation:** Add compile-time assertions for `std::ranges::forward_range`, iterator copyability and non-writability. Add negative assertions for unsupported borrowed, contiguous, random-access or bidirectional concepts. Runtime-test copied iterator independence, repeated traversal and nested traversal. Destroy, move and assign the owning range, then verify that previously extracted iterators are rejected or remain outside the supported contract as documented.

### AC-0036-03 — Whole-store enumeration

- **Requirement:** `AllEntities<K>(shape)` returns every live entity exactly once, including entities in disconnected components and entities not reachable from an enumerated root.
- **Complexity:** Document and measure O(S) for sparse-slot scanning unless a maintained live index provides a stronger bound. No query may silently claim O(L).
- **Owner:** Topology store and range implementation owners.
- **Validation:** Test dense, sparse, all-deleted and empty stores. Verify no tombstones, duplicate IDs or cross-snapshot IDs. Record construction and iteration time separately.

### AC-0036-04 — Local traversal

- **Requirement:** `CoedgesOf` and `RadialCoedgesOf` preserve their anchor, relation direction, multiplicity and finite closure semantics.
- **Complexity:** For k emitted uses, total traversal is O(k), with O(1) traversal state after snapshot retention.
- **Allocation:** Zero allocations and zero ownership-reference increments per step.
- **Owner:** Topology traversal implementation owner.
- **Validation:** Test one-item, long, sparse, seam and malformed fixtures. Instrument every increment and dereference. Pass only if measured per-step allocation and ownership traffic are zero.

### AC-0036-05 — Distinct adjacency preparation

- **Requirement:** `FacesOf` and `EdgesOf` never hide a whole-store scan, hash-set allocation or preparation object inside ordinary local range construction.
- **Complexity:** Preparation and memory are explicitly documented and measured. Query iteration is deterministic and uses ascending entity index after preparation/materialisation.
- **Owner:** Topology query-preparation owner.
- **Validation:** Compare maintained-index, explicit-preparation and eager-result paths. Instrument scans, allocations and preparation-token retention. Reject any implementation whose local query cost depends on unrelated model size without an explicit contract.

### AC-0036-06 — Recursive result

- **Requirement:** `CollectDescendants<K>` returns a complete `OwnedTraversalResult<K>` or `KernelError`; it never returns a partial successful result.
- **Complexity:** Unique traversal targets O(V+A) time and O(V) visited state, excluding documented sorting. Output-sensitive occurrence traversal documents its output bound.
- **Owner:** Recursive traversal implementation owner.
- **Validation:** Test root inclusion, exclusion, deduplication, disconnected roots, cycles and invalid category preconditions. Verify deterministic depth-first order and complete-result error behavior.

### AC-0036-07 — Prepared-index lifetime

- **Requirement:** `TopologyQueryPreparation` retains the snapshot token and remains immutable after construction. Prepared data remains valid while any preparation or result references it and is reclaimable afterward.
- **Complexity:** Preparation ownership work is O(1) per retained snapshot token, excluding documented index construction.
- **Owner:** Topology query-preparation owner.
- **Validation:** Use reference-count instrumentation and weak lifetime probes. Create, copy, move and destroy preparations and results in every order. Verify exactly one logical snapshot remains alive while any owner exists and that optional prepared indexes are reclaimed after the final owner releases them.

### AC-0036-08 — Benchmark and review gate

- **Requirement:** Measurements cover construction, iteration, allocations, ownership traffic and retained memory.
- **Owner:** Topology performance owner.
- **Validation:** Run dense and sparse stores, high-valence vertices, long radial rings, multiple page counts and recursive workloads. Attach results to the ADR review. Status remains Proposed until all acceptance criteria pass and integration review confirms compliance with snapshot, topology and preparation ownership decisions.

## Consequences

Algorithms gain a reusable, range-compatible interface and cannot accidentally mutate connectivity while iterating. We must implement and test distinct query families, incidence preparation and ownership. Allocation-free local traversal does not imply allocation-free deduplication or recursive search.

C/Python adapters may wrap these contracts using opaque cursors or materialised values under ADR-0020; they must preserve multiplicity and lifetime semantics. This ADR does not add STL types to the C ABI or freeze binding signatures.


## Normalization and remaining implementation work

The following conventions are accepted across ADR-0036 through ADR-0038:

- **Whole-store entity query:** `AllEntities<K>(shape)`.
- **Local entity queries:** `LoopsOf`, `FacesOf` and `EdgesOf`.
- **Use queries:** `CoedgesOf` and `RadialCoedgesOf`; these preserve incidence multiplicity.
- **Recursive result:** `CollectDescendants<K>(shape, root, options)` returns `std::expected<OwnedTraversalResult<K>, KernelError>`.
- **Fallible acquisition:** Every public query acquisition uses `std::expected<Range, KernelError>`.
- **Diagnostic validation:** Validation failures are reported as structured `KernelError` values, never as ordinary end-of-range.
- **Diagnostic cursor outcomes:** `Item`, `End` and `Error`.
- **Stored connectivity identity:** `CoedgeId`.
- **Oriented traversal value:** `OrientedCoedgeUse`.
- **Orientation terminology:** `stored edge-use sense`, `traversal reversal` and `effective edge direction`.
- **Manifold helper:** `ManifoldOppositeFace`, with explicit precondition validation.
- **Snapshot terminology:** `snapshot`, `snapshot ownership token` and `snapshot descriptor`.
- **Mutation invalidation terminology:** `topology epoch`; entity generations do not validate builder cursors.
- **Builder ownership terminology:** stable `BuilderToken` plus topology epoch.
- **Prepared-index ownership:** maintained indexes belong to the snapshot; explicit prepared indexes belong to `TopologyQueryPreparation`.
- **Result terminology:** `Range` for lazy forward ranges, `OwnedTraversalResult<K>` for eager recursive results and `KernelError` for fallible operations.
- **Recursive relation order:** solid-to-shell, shell-to-face, face-to-loop, loop-to-coedge, coedge-to-edge and edge-to-vertex.

The following items are implementation details rather than unresolved ownership or validation decisions:

1. The concrete C/Python adapter spellings and opaque-handle layout.
2. The separate geometric-sector query API.
3. The future general occurrence-path type.
4. Benchmark thresholds for range construction, local traversal, preparation, retained memory and ownership traffic.
5. The exact internal reference-counting primitive used by the snapshot control block, provided it satisfies the ownership and complexity contracts above.

## Integration acceptance criteria

### AC-INT-01 — ADR allocation and repository integration

- **Requirement:** Recheck ADR numbering, split this document into the three proposed files and add entries to `docs/adr/README.md`.
- **Owner:** ADR maintainer.
- **Validation:** Confirm no numbering collision, correct filenames, valid index links and preserved `Proposed` status. CI must verify that every referenced ADR file exists and every index entry resolves.

### AC-INT-02 — Reciprocal references

- **Requirement:** Add reciprocal references in ADR-0004, ADR-0005, ADR-0018 and ADR-0027 for adjacency views, snapshot token lifetime, independent bounded validation and cursor epoch invalidation.
- **Owner:** ADR maintainer.
- **Validation:** Review each target ADR for a reciprocal link and a concise statement of the dependency. Existing status fields must remain unchanged.

### AC-INT-03 — Foundational implementation order

- **Requirement:** Implement `TopologySnapshot`, `ValidationCertificate`, `BuilderToken` and topology epochs before public ranges.
- **Owner:** Topology architecture owner.
- **Validation:** Enforce dependency order in the implementation plan and code review. Public range code must not compile against a borrowed Shape-only lifetime model or uncertified topology.

### AC-INT-04 — Staged feature allocation

- **Requirement:** Implement local entity, loop and radial ranges with snapshot retention first; implement maintained incidence indexes and explicit distinct-adjacency preparation next; implement recursive discovery only when a concrete consumer is identified.
- **Owner:** Topology implementation lead.
- **Validation:** Track each stage as a separately reviewed milestone with passing ADR-specific acceptance criteria. No recursive API may be marked complete without a named consumer, complexity measurement and complete-result error tests.

### AC-INT-05 — Ownership and validation contract preservation

- **Requirement:** Existing Arena behavior, implementation convenience and benchmark results must not silently replace the resolved ownership or validation contracts.
- **Owner:** Topology architecture owner and reviewer.
- **Validation:** Require explicit ADR amendment for any deviation. Code review must reject page-table copying per range, hidden whole-store scans, builder commit as certification, implicit preparation ownership or corruption encoded as end-of-range.

### AC-INT-06 — Proposed-status exit gate

- **Requirement:** ADR-0036, ADR-0037 and ADR-0038 remain Proposed until API sketches, ownership implementation, validation diagnostics, iterator-lifetime tests, integration references and benchmark results have been reviewed.
- **Owner:** ADR approver.
- **Validation:** Record a signed review checklist containing all AC-0036, AC-0037, AC-0038 and AC-INT criteria. Change status only after every mandatory criterion passes or an explicit follow-up ADR is accepted.

## Inspected sources

All links are pinned to the reviewed commit:

- [ADR template](https://github.com/onurtuncer/Stereon/blob/feb97bf4222be24d3919141245b47510fc7aec97/docs/adr/0000-template.md)
- [Topology storage and handles](https://github.com/onurtuncer/Stereon/blob/feb97bf4222be24d3919141245b47510fc7aec97/docs/adr/0004-topology-storage-and-handles.md)
- [Immutability and sharing](https://github.com/onurtuncer/Stereon/blob/feb97bf4222be24d3919141245b47510fc7aec97/docs/adr/0005-immutability-and-sharing.md)
- [Error handling](https://github.com/onurtuncer/Stereon/blob/feb97bf4222be24d3919141245b47510fc7aec97/docs/adr/0006-error-handling.md)
- [Threading and determinism](https://github.com/onurtuncer/Stereon/blob/feb97bf4222be24d3919141245b47510fc7aec97/docs/adr/0007-threading-and-determinism.md)
- [Naming conventions](https://github.com/onurtuncer/Stereon/blob/feb97bf4222be24d3919141245b47510fc7aec97/docs/adr/0017-naming-conventions.md)
- [Controlled topology mutation](https://github.com/onurtuncer/Stereon/blob/feb97bf4222be24d3919141245b47510fc7aec97/docs/adr/0027-controlled-topology-mutation-and-euler-operators.md)
- [Current topology documentation](https://github.com/onurtuncer/Stereon/blob/feb97bf4222be24d3919141245b47510fc7aec97/docs/topology.rst)
- [Current typed-handle implementation](https://github.com/onurtuncer/Stereon/blob/feb97bf4222be24d3919141245b47510fc7aec97/src/core/Handle.cppm)
