# ADR-0027: Controlled Topology Mutation and Euler Operators

- **Status:** Proposed
- **Date:** 2026-10-13
- **Owner:** Onur Tuncer
- **Phase:** 0 (contract); implementation staged with topology work
- **Related:** ADR-0004, ADR-0005, ADR-0006, ADR-0009, ADR-0011, ADR-0016, ADR-0018, ADR-0025, ADR-0026

## Context

Stereon's topology store uses typed handles, ordered coedge loops, radial rings, and mutable builders that produce immutable shapes (ADR-0004 and ADR-0005). These structures contain reciprocal relationships: changing one link usually requires coordinated changes to other links, ownership, geometry associations, and derived data.

If modelling algorithms update those records independently, each algorithm must implement the same connectivity rules. A locally plausible edit can leave a broken loop, an incomplete radial ring, a stale primary coedge, or an invalid geometric association. Such faults become harder to diagnose after several operations have composed.

Euler operators provide a useful basis for controlled topological construction and modification. Their applicability depends on the supported model category and their stated assumptions. An Euler characteristic alone does not prove that a representation is a valid manifold or a geometrically valid solid, and classical manifold operators do not establish a complete operator set for Stereon's non-manifold intermediates.

Stroud's *Boundary Representation Modelling Techniques* motivates examining Euler operators and stepwise modelling. The policy below is Stereon's own design decision; it does not adopt an unverified catalogue or algorithm from the book.

## Options

1. **Allow each algorithm to edit topology records directly.** Minimal initial infrastructure, but connectivity logic and failure handling are duplicated across algorithms.
2. **Require classical Euler operators for every edit.** A disciplined basis for applicable model categories, but insufficient as an unqualified rule for radial non-manifold structures, degeneracies, and geometry-association changes.
3. **Use a controlled topology mutation layer with contracted operators.** Euler operators are included where their assumptions apply; other edits have equally explicit contracts. Private transactions compose primitives into valid public results.

## Decision

Adopt option 3.

All production changes to topology connectivity shall pass through a controlled mutation layer on `ShapeBuilder`. Its operators own the corresponding coordinated record changes. Modelling algorithms shall not bypass this layer by writing incidence, ownership, loop, or radial-ring fields directly.

This is a semantic and access-control boundary, not a requirement to expose low-level operators through the public modelling API.

### 1. Mutations have explicit contracts

Each operator shall specify:

- Accepted entity kinds and topology categories, including manifold, sheet, wire, and non-manifold restrictions.
- Required handle/store validity and local connectivity conditions.
- Which entities and incidence relations may change.
- The structural invariants preserved and any invariants deliberately suspended inside a private transaction.
- Required geometric-association updates and invalidation of certificates or derived data.
- Output handles, creation/deletion/modification records, and failure behavior.
- The applicable local checks and requirements before transaction commit.

External data errors and unsupported configurations produce structured errors under ADR-0006. Internal handle misuse and violated programmer obligations remain contracts under ADR-0016. Operators must not turn invalid imported geometry into a host-terminating contract violation.

An operator's name is not its specification. Names and signatures are fixed when the first implementation and its contract are reviewed.

### 2. Separate connectivity edits from geometric construction

Topology operators manipulate connectivity; geometric algorithms compute points, curves, surfaces, pcurves, parameter intervals, and certificates. The association boundary follows ADR-0025.

Some complete edits require both. Splitting an edge, for example, requires a consistent split across its incident coedges, parameter restrictions, geometric representations, and endpoint associations. A connectivity-only edit may leave these obligations pending inside a private transaction, but it cannot publish a completed B-rep result.

ADR-0009 remains authoritative: primary coedges are explicit; radial-ring reordering must not change edge evaluation; surface-defined curves are not replaced silently by independent 3D curves. An edit that invalidates certificates must regenerate them or fail before committing a result that requires them.

### 3. Euler operators are a contracted subset

Use Euler operators for construction and modification where their mathematical preconditions match the represented topology. Document each operator's effects on entity counts, connectivity, connected components, and genus where those quantities apply.

Do not use a single global formula such as `V - E + F = 2` as a universal acceptance criterion. Faces with multiple loops, multiple shells, cavities, disconnected components, wires, and non-manifold structures require the appropriate representation-specific accounting and checks.

The first operator catalogue shall identify its supported domain and demonstrate representative construction sequences. This ADR does not claim completeness or require that every admissible B-rep can already be generated.

Every paired inverse shall state the conditions under which it is an inverse. Restoring topology up to an explicit entity correspondence does not require reusing deleted IDs or reproducing the original storage layout.

### 4. Radial and ownership edits are controlled

The mutation layer shall provide the operations needed to maintain reciprocal loop links, edge-use radial rings, and parent ownership together. It must support the non-manifold representation of ADR-0004 rather than hard-code two incident face uses for every edge.

For each edit, affected loops and radial rings must have unambiguous membership. A coedge cannot silently belong to two loops or refer to an unrelated edge in its radial ring.

Deleting or replacing a primary coedge requires an explicit reassignment consistent with ADR-0009's evaluation semantics, with affected associations revalidated. Reordering a ring alone is not permission to select a new primary.

These are invariants, not a prescribed list of setter functions. Unrestricted public setters for individual connectivity fields do not satisfy this decision.

### 5. Private transactions bound temporary inconsistency

A compound topology edit executes within a private builder transaction. Before it commits, intermediate records are inaccessible to ordinary model queries, serialization, meshing, and publication as a valid `Shape`.

Internal algorithms may inspect intermediate states only through interfaces whose contracts admit those states. Pending geometry obligations and suspended structural invariants must be tracked explicitly; an undocumented half-edited graph is not a supported intermediate representation.

On recoverable failure, either:

- the transaction restores the builder's entry state; or
- the private candidate builder is discarded and cannot be frozen or reused as a valid result.

A reusable caller-visible builder must roll back if its API promises continued use. Discard is appropriate for a private candidate owned by one modelling operation. The implementation must choose and document the behavior rather than silently return a damaged reusable builder.

These transactions stage private work; they do not introduce a public shape undo log. Frozen inputs and their shared pages remain unchanged on success and failure under ADR-0005 and ADR-0026.

### 6. Commit and publication have explicit gates

Local post-checks run at operator or transaction boundaries as appropriate. They cover the touched relations and their necessary incidence neighborhood, not just the records directly written.

A transaction may commit only when all invariants promised by its contract are restored. Publication of an operation result additionally requires the checks and guarantees of ADR-0018 and ADR-0026. In particular:

- No unfinished connectivity edit or unresolved association obligation may escape.
- The applicable structural, geometric, and solid-level checks follow ADR-0018's build and Context policies.
- Certificates required by ADR-0009 must be present and current.
- Shape, history, tolerance changes, and diagnostics must be published consistently in `OperationResult`.

Local checks supplement the independent validity checker; they do not replace it or reuse its independence requirement as a reason to duplicate all checks on every primitive edit.

### 7. Identity, history, and ordering follow existing decisions

Creation, deletion, and modification events shall be available to the operation's history recorder (ADR-0011). Low-level events are not automatically a complete semantic history; the enclosing modelling operation remains responsible for meaningful source-to-result relations.

Do not reuse slots deleted in the same builder session. `Freeze()` preserves indices and tombstones; `Compact()` remains the explicit renumbering operation, as specified by ADR-0004. Shape-local handles do not establish persistent identity across operations.

Operator composition and history output shall respect ADR-0007's deterministic ordering. Private topology transactions do not make a builder concurrently mutable; ADR-0004's single-threaded builder policy remains in force.

### 8. Loading and repair do not create hidden bypasses

Import and repair operations must uphold the same published-state guarantees. A bulk loader may populate private tables through a dedicated controlled construction path with a declared validation and commit contract; it must not expose unchecked tables as valid topology.

Tests may deliberately construct invalid records through isolated fault-injection facilities. Such facilities are unavailable to production modelling paths and exist to test the checker and failure handling.

## Consequences

- Connectivity rules and coordinated edits have identifiable owners and reviewable contracts.
- Complex modelling operations can compose smaller verified mutations.
- Euler operators provide useful construction and inverse-test cases without imposing manifold-only assumptions on the entire store.
- Transactions and association tracking add implementation work and may increase temporary storage.
- Local validation requires careful neighborhood tracking to avoid both missed dependencies and repeated whole-model scans.
- The independent validity checker remains necessary because the mutation layer can itself contain bugs.
- The operator catalogue can grow incrementally; unsupported edits are reported explicitly until implemented.

## Verification

The first implementation slice shall establish the mutation boundary and demonstrate:

1. **Controlled access:** production algorithms cannot directly write protected connectivity fields; internal storage access is confined to the mutation and controlled-loading implementations.
2. **Representative edits:** vertex/edge construction, edge splitting, loop/face splitting where supported, and radial membership changes preserve their declared local invariants.
3. **Domain coverage:** valid manifold, open-sheet, wire, and at least one radial non-manifold configuration exercise the relevant supported primitives; unsupported cases fail explicitly.
4. **Inverse properties:** valid operator/inverse pairs restore connectivity up to explicit correspondence and applicable geometric agreement, without assuming recycled handles.
5. **Failure containment:** injected failures after meaningful stages either restore the reusable builder or discard the private candidate. Previously frozen input shapes remain unchanged through all aliases.
6. **Association consistency:** edge splits update every affected coedge and parameter restriction; changing a pcurve invalidates dependent certificates; no uncertified required association is published.
7. **Primary-side stability:** radial reordering leaves evaluation unchanged, while primary-coedge deletion either performs a validated explicit reassignment or fails.
8. **Independent checking:** completed examples pass the applicable ADR-0018 levels; seeded broken links, ownership errors, and incomplete rings are rejected.
9. **Composition and determinism:** a multi-step edit that fails late publishes nothing; successful public operations satisfy ADR-0007 G1.
10. **Property-based coverage:** generated valid operator sequences satisfy local invariants and final independent validation, as required for topology code by ADR-0015.

Tests are implemented with the corresponding operators. These are acceptance requirements, not a claim that an operator implementation or complete catalogue already exists.

## References

- Ian Stroud, *Boundary Representation Modelling Techniques*, Springer, 2006, [DOI: 10.1007/978-1-84628-616-2](https://doi.org/10.1007/978-1-84628-616-2). Relevant chapters: “Euler operators” and “Stepwise modelling algorithms.” Publisher-level coverage has been verified; detailed chapter review is pending. No specific operator algorithm is attributed to an unread passage.
- [ADR-0004: Topology storage and handles](0004-topology-storage-and-handles.md)
- [ADR-0005: Immutability and structural sharing](0005-immutability-and-sharing.md)
- [ADR-0006: Error handling](0006-error-handling.md)
- [ADR-0009: Edge curves as surface-intersection curves](0009-parameter-space-representation.md)
- [ADR-0016: C++26 adoption and contracts policy](0016-cpp26-adoption-and-contracts.md)
- [ADR-0018: Validity checking](0018-validity-checking.md)
- [ADR-0025: Separate Topology from Geometry](0025-separate-topology-from-geometry.md)
- [ADR-0026: Specify Modelling Operations by Contracts](0026-modelling-operations-by-contracts.md)

This proposal supplements the existing decisions and supersedes none. It specifies a controlled mutation boundary without selecting a complete operator catalogue, changing the storage representation, or accepting the related proposals.
