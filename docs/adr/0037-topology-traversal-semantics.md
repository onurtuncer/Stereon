# ADR-0037: Entity, use and orientation semantics for topology traversal

- **Status:** Proposed
- **Date:** 2026-10-17
- **Owner:** Onur Tuncer
- **Phase:** Topology foundation
- **Related:** [ADR-0004](0004-topology-storage-and-handles.md), [ADR-0007](0007-threading-and-determinism.md), [ADR-0009](0009-parameter-space-representation.md), [ADR-0011](0011-history-and-naming.md), [ADR-0018](0018-validity-checking.md), [ADR-0025](0025-separate-topology-from-geometry.md), [ADR-0027](0027-controlled-topology-mutation-and-euler-operators.md), [ADR-0030](0030-entity-lineage-in-operation-results.md), [ADR-0036](0036-typed-topology-ranges.md), [ADR-0038](0038-topology-traversal-lifetime-and-validation.md)

## Context

An edge identity is not its use in a boundary. A periodic face may use one seam edge twice, with different coedges and pcurves. A non-manifold edge may have more than two incident face uses. Removing duplicates or assuming a unique opposite face can therefore erase required topology.

ADR-0004 supplies ordered loop links and radial coedge rings. Their physical order alone does not define geometric angular order or occurrence orientation through a parent hierarchy.

## Options

1. **Deduplicate all results automatically.** Convenient entity counting; destroys boundary multiplicity.
2. **Return only raw uses.** Preserves topology; forces every identity query to implement deduplication.
3. **Separate entity queries from use queries and define order/orientation explicitly.** More names, but unambiguous algorithm inputs.

## Decision

Adopt option 3.

### Entity versus use identity

An entity query returns each relevant shape-local entity once. A use query returns each represented incidence, even when several uses reference the same underlying entity. There is no implicit conversion from an occurrence traversal to a deduplicated traversal.

For edge uses, retain at least the coedge, underlying edge, owning loop/face context and traversal sense needed to recover the correct pcurve and parameter direction. A use object references existing data; it does not duplicate geometry or invent a new persistent identity. For repeated higher-level uses, preserve parent-use/path context when the logical model permits repetition. A general occurrence representation is not established by this ADR.

`FacesOf(shape, edge)` deduplicates by FaceId within that snapshot. `RadialCoedgesOf(shape, edge)` does not deduplicate by face, edge or geometric coincidence. Geometrically coincident entities remain distinct unless a modelling operation explicitly identifies them.

### Finite loop traversal

A closed loop is exposed as a finite range: start at its designated anchor, follow Next, emit each coedge once and stop before emitting the anchor again. A one-coedge loop emits one item. An empty loop may be returned only where that topology category permits it; an illegal empty boundary is a validation error.

An explicitly reversed walk follows Prev and reverses each emitted use’s traversal sense. With the same anchor c0, a stored cycle `[c0,c1,c2]` becomes `[reverse(c0),reverse(c2),reverse(c1)]`. This is anchored reversal, not necessarily the exact sequence produced by a generic container reverse adaptor. Double reversal restores the original sequence and senses.

A plain CoedgeId identifies stored connectivity. An OrientedCoedgeUse adds traversal sense; it does not replace the record’s relationship to its underlying edge. `StoredEdgeUseSense` returns the authoritative stored sense. `EffectiveEdgeDirection` composes that sense with `TraversalReversal` exactly once. Never apply the same reversal twice or manually mask packed bits.

When an API accepts an oriented parent use, its contract must name the orientation reference frame and compose each relevant incidence sense once. Raw entity enumeration has no inherited orientation. Geometric normal evaluation and assembly placement are outside this decision.

### Radial traversal and adjacency

A radial walk emits all coedges in the stored ring exactly once, starting from a designated radial anchor. That anchor need not be the edge’s evaluation primary coedge. Reordering or reanchoring the ring must not change the primary coedge or edge evaluation semantics of ADR-0009.

Ring order is a deterministic connectivity order for a fixed snapshot. It is **not** a promise of angular ordering around a tangent. Algorithms needing geometric sector order must request a separate geometric query, which is outside this ADR and must define degeneracy and error handling before implementation.

Do not provide an unqualified `Twin` or `OppositeFace` operation for arbitrary topology. `ManifoldOppositeFace` validates its manifold preconditions and returns an explicit error when they fail; ring size two alone is not sufficient proof of the desired manifold neighbourhood.

|Case                                           |Required behaviour                                                                   |
|-----------------------------------------------|-------------------------------------------------------------------------------------|
|Seam edge used twice by one face               |Two coedge uses, correct pcurve for each; one distinct face                          |
|Edge shared by three faces                     |All three uses and all three distinct faces; no forced pairing                       |
|Boundary edge with one use                     |One radial item; no invented opposite face                                           |
|Spatial/wire edge with no face use             |Empty radial/face query if permitted by the edge category                            |
|Closed edge with coincident endpoint identities|Unique vertex/edge adjacency deduplicates; endpoint-use queries preserve both roles  |
|Degenerate edge                                |Preserve the topological entity and permitted uses; do not filter by geometric length|

### Ordering and identity

- Whole-store entity ranges use ascending live slot index.
- Loop/radial ranges preserve their documented anchor and links.
- Distinct adjacency results use ascending entity index after explicit preparation/materialisation; hash iteration order is never observable.
- Recursive entity collection uses deterministic depth-first discovery, first-visit deduplication and the fixed relation order defined in this ADR: solid-to-shell, shell-to-face, face-to-loop, loop-to-coedge, coedge-to-edge and edge-to-vertex.
- Radial incidence, reciprocal links and vertex incidence are not traversed as containment relations.
- No guarantee equates traversal order across `Compact()`, separate imports, different construction histories or geometrically equivalent models. Comparisons across outputs require remaps/history. No iterator position, slot index or occurrence path becomes a persistent name under ADR-0011/0030.

## Resolved supporting decisions

### Recursive ordering decision

`CollectDescendants<K>` uses deterministic depth-first discovery with first-visit deduplication. The fixed relation order is:

1. Solid to shells;
2. Shell to faces;
3. Face to loops;
4. Loop to coedges;
5. Coedge to edges;
6. Edge to vertices.

Relations not listed above are not treated as containment relations. In particular, radial incidence, reciprocal links and vertex incidence are not traversed as descendants. Within each relation, the query uses the local ordering contract of the source range. Root inclusion, deduplication and ordering are controlled by `CollectDescendantsOptions`; the default includes the root when its kind matches K, deduplicates by entity identity and uses the ordering above.

### Occurrence representation decision

No general persistent occurrence-path type is introduced. `OrientedCoedgeUse` may retain coedge, edge, parent loop/face context and traversal sense. Recursive results may optionally retain transient parent-use context when required by a consumer, but that context is not an entity identity and is not serialized as a stable name. A future occurrence-path ADR is required before exposing occurrence paths as a public cross-binding type.

## Acceptance criteria

### AC-0037-01 — Entity/use multiplicity

- **Requirement:** Entity queries return each relevant entity once; use queries preserve every represented incidence.
- **Owner:** Topology semantics implementation owner.
- **Validation:** Use fixtures containing seam reuse, non-manifold edges, boundary edges, spatial edges, closed edges and degenerate edges. Compare exact expected entity IDs, use counts and parent contexts. Pass only if no use is removed by implicit deduplication and no coincident entities are merged.

### AC-0037-02 — Loop order and reversal

- **Requirement:** Forward loop traversal emits the anchored stored cycle exactly once. Reversed traversal follows Prev from the same anchor and reverses each emitted use exactly once.
- **Complexity:** O(k) for k emitted coedges, with O(1) traversal state.
- **Owner:** Loop traversal implementation owner.
- **Validation:** Test triangle, one-coedge, hole, seam and reversed fixtures. Verify exact sequences, directed endpoints, continuous orientation where required and double-reversal identity. Include a test proving that generic container reversal is not substituted for anchored reversal.

### AC-0037-03 — Radial order and primary coedge

- **Requirement:** Radial traversal emits every coedge in the stored ring exactly once and does not redefine the edge’s primary coedge or evaluation semantics.
- **Owner:** Radial traversal and edge-evaluation owners.
- **Validation:** Reanchor and reorder radial rings in equivalent fixtures. Verify exact ring output for each snapshot and unchanged primary-coedge selection and evaluated edge geometry. Test one-, two- and three-or-more-use rings.

### AC-0037-04 — Orientation composition

- **Requirement:** `EffectiveEdgeDirection` composes stored edge-use sense and traversal reversal exactly once.
- **Owner:** Orientation API owner.
- **Validation:** Exhaustively test all combinations of stored sense and traversal reversal, including nested oriented parent uses. Verify no direct packed-bit access in the public implementation through static analysis or code review. Confirm raw entity enumeration does not inherit orientation.

### AC-0037-05 — Manifold helper preconditions

- **Requirement:** `ManifoldOppositeFace` succeeds only for exactly two valid coedge uses belonging to distinct faces and satisfying all declared manifold invariants.
- **Owner:** Manifold-query owner.
- **Validation:** Test boundary, seam, non-manifold, duplicate-face and invalid-ring cases. Verify explicit `KernelError` values and absence of unqualified `Twin` or `OppositeFace` APIs. Ring size two alone must not make the query succeed.

### AC-0037-06 — Deterministic ordering

- **Requirement:** Whole-store, local, distinct-adjacency and recursive queries produce the documented deterministic order for a fixed snapshot.
- **Owner:** Topology query implementation owner.
- **Validation:** Repeat queries, vary scheduler execution order and run under thread sanitizers. Compare exact sequences. Verify distinct results are ascending by entity index and recursive results follow fixed depth-first relation order. Compare compacted snapshots only through explicit remaps.

### AC-0037-07 — Acceptance fixture coverage

- **Requirement:** The fixture suite includes a triangle loop, one-coedge loop, face with a hole, periodic seam, open sheet, disconnected shells, spatial edge, closed edge, degenerate edge and an edge with at least three face uses.
- **Owner:** Topology test owner.
- **Validation:** Store independently recorded expected sequences, multiplicities, identities and senses. The ADR remains Proposed until every fixture passes all applicable criteria.

## Consequences

Boundary algorithms retain seams and orientation, while analysis/counting queries can request unique entities. Non-manifold topology is represented without inventing manifold neighbours. Clients must choose the query matching their mathematical task. Canonical geometric ordering and persistent naming remain separate problems.

