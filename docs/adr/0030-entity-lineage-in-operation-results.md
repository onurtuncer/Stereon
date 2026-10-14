# ADR-0030: Make Entity Lineage Part of Operation Results

- **Status:** Proposed
- **Date:** 2026-10-14
- **Owner:** Onur Tuncer
- **Phase:** 3 (v0.3 construction); composition through v0.4–v0.6 Booleans
- **Extends:** ADR-0004, ADR-0006, ADR-0011, ADR-0026, ADR-0027

## Context

ADR-0011 selects operation history over full persistent naming. The four relation labels alone do not settle how to qualify input references, describe one-to-many splits and many-to-one merges, distinguish a generator from a surviving entity, or represent primitives that have no input topology.

Storage indices cannot fill these gaps. Two inputs can both contain `FaceId(7)`, shared immutable geometry may belong to different topology, and compaction can renumber an otherwise unchanged entity. Downstream attributes and shape sensitivities need semantic correspondence with explicit limitations.

Stroud §8.3 discusses recording changed entities and preserving information during split/merge operations. The relation model and publication guarantees below are Stereon's design, not a persistent-naming solution claimed by the book.

## Options

1. **Infer identity from equal indices, addresses or geometry.** Cheap but confuses independent inputs, duplication, coincidence and compaction.
2. **Record only low-level creation/deletion events.** Useful debugging data; loses the semantic relationship of a split, merge or generated face.
3. **Publish an operation-scoped lineage relation with qualified references and explicit event semantics.** Additional storage and recorder work; supports reliable downstream interpretation.

## Decision

Adopt option 3 as the detailed contract for ADR-0011's `History` in `OperationResult`. History is a many-to-many relation, not a function from each input to exactly one output. It describes one operation on particular immutable snapshots; it does not establish identity across independent regeneration.

### 1. References are qualified by their owning snapshot

Every source reference identifies the operation input ordinal, entity kind and shape-local index; every target identifies the output snapshot/result ordinal, kind and index. These ordinals refer to an immutable input/output manifest held by the result or explicitly retained history resource. Bare index equality is never sufficient.

The manifest distinguishes multiple argument occurrences even when they share storage. Operation-local references are qualified by the owning result/history; they are not process-global IDs. Standalone persisted histories carry serializable snapshot/provenance references rather than raw addresses. Native load must not silently bind a reference to the current shape because its index happens to match.

Accessors document ownership. Keeping an `OperationResult` or retained `History` keeps the manifests and required snapshots alive for dereferenceable references; a separately exported provenance-only record may contain detached references but must report that they cannot be dereferenced until their snapshots are supplied. Avoid accidental lifetime chains beyond the snapshots explicitly retained by the API.

Bindings follow this rule too. ADR-0020's `{kind, index}` entity reference remains relative to a shape resource; history accessors additionally expose its input/output owner. Do not change the compact topology handle layout to encode global identity.

### 2. Separate continuation, generation and deletion

| Relation / event | Meaning |
| --- | --- |
| `Unchanged` continuation | The input entity survives with the same promised geometric/topological meaning; a storage remap alone does not make it modified. |
| `Modified` continuation | An output is a semantic continuation or retained portion of an input, with a documented change. Includes split pieces and merged continuations. |
| `Generated` provenance | An entity is newly created by the operation from named input generators or an operation root. Being a generator does not imply survival. |
| `Deleted` disposition | The input has no surviving continuation in the returned model. It may still have generated other entities. It has no fake/null output entity. |
| Split event | One source has several retained continuations identified explicitly by the operation, with output members listed. |
| Merge event | Several sources have one retained continuation identified explicitly by the operation, with source members listed. |

Keep ADR-0011's four relation labels. Split and merge are explicit grouped semantic events over continuation relations, not replacement labels and not guesses from the count of all provenance edges. An edge extruded into a face has generated a face; it has not necessarily split into one. Complex many-to-many changes can contain both split and merge events; record the known relation without forcing a bijection.

Generators may have a different entity kind from their outputs. A fillet face can name an edge and supporting faces. Primitive construction and import with no input entities use a creation root naming the operation/parameters or source document and record identity where available. They do not fabricate an input entity to satisfy completeness.

`Unchanged` is not inferred solely from a retained slot or shared geometry. Any metadata-only changes are reported through the metadata diagnostics/policy contract even when the underlying geometry/topology continues unchanged.

### 3. History is complete for the declared participating domain

Every returned live entity has an explicit continuation, generation provenance or a creation root. Every participating live input entity has a continuation or a Deleted disposition; absence from a sparse table must not ambiguously mean either deletion or unchanged.

An implementation may encode unchanged ranges or shared substructures compactly, provided the public history query is complete. An operation on a selected subshape states its participating domain; entities outside it are not falsely marked deleted. A solid Boolean's input domain normally includes both complete operands, including the tool.

A legitimate empty result still carries history: all participating entities lacking surviving continuations have Deleted dispositions. If a removed face generates an intersection edge, both deletion and generation remain visible.

The enclosing operation owns these semantics. Low-level topology-mutation events from ADR-0027 are evidence for the recorder, not automatically the final public history. Temporary intermediate entities that are created and removed privately are not published as live result entities.

### 4. Compose relations without inventing identity

Composition matches intermediate references only when their immutable snapshot manifests agree. Reject incompatible manifests; equal local indices in unrelated results are insufficient. Preserve generator ancestry separately from continuation paths and preserve deletion semantics when no final continuation remains.

Retain event provenance needed to explain split-then-merge and other many-to-many sequences. A composed lineage view can summarize end-to-end descendants, but must not claim the middle events never occurred or turn every shared ancestor into an identity relation. Record deterministic ordering by input ordinal, kind, index and stable operation-local event key; never by pointer value or task completion order.

`Freeze()` keeps indices/tombstones under ADR-0004. `Compact()` produces a total remap for live entities, with semantic Unchanged continuations despite different indices; compose that remap into references. Builder generations protect stale builder handles and are not persisted naming IDs.

### 5. Publish history atomically with the model

A successful `OperationResult` contains mutually consistent shape, history, tolerance changes and metadata diagnostics. All published targets refer to live entities of the returned snapshot. Missing mandatory history or invalid references fail the operation before publication, following ADR-0026.

On failure, immutable inputs and already published results remain unchanged. Diagnostic traces may describe attempted private edits, but must not appear as committed history. Attribute propagation consumes this validated semantic lineage under ADR-0031.

A lineage relation alone supplies neither a geometric point/parameter correspondence nor derivatives. Operations offering shape sensitivities must additionally supply the correspondence and regularity information required by ADR-0021. Persistent naming across regeneration remains outside the v1.0 promise.

## Consequences

Lineage becomes usable across multiple inputs and chained operations. Retained snapshots and event provenance consume memory; compact unchanged encodings and explicit resource lifetimes limit the cost. Algorithms must distinguish semantic continuation from construction dependencies, and history needs its own consistency checks rather than relying on shape validity alone.

## Verification

1. Two inputs with equal local face indices remain distinct in history, including two argument occurrences sharing storage; references survive release of the caller's original shape resources when retained by the result.
2. Splitting one edge/face reports every retained continuation and its split group. Merging multiple sources reports every contributor and its merge group; generated side faces are not mislabelled as splits.
3. Primitive construction and import without input topology report creation roots. A deleted generator can still have Generated outputs. Empty successful results expose complete deletion dispositions.
4. Split → modify → merge → compact chains compose to the correct endpoint lineage while retaining intermediate event provenance. Unrelated intermediate manifests are rejected.
5. Compaction changes indices but not geometric/topological continuation semantics. No returned target points at a tombstone, discarded builder or private temporary entity.
6. Inject a failure during history assembly or final publication: no partial model/history escapes. Run G1 checks on canonical histories as well as geometry across thread counts.
7. Attribute and derivative consumers demonstrate that they distinguish provenance from same-kind continuation and refuse missing correspondence rather than inventing it.

## References and scope

- Ian Stroud, *Boundary Representation Modelling Techniques*, Springer, 2006, §8.3, printed pp. 262–270; information handling and split/merge discussion reviewed in the supplied book.
- [ADR-0004](0004-topology-storage-and-handles.md), [ADR-0006](0006-error-handling.md), [ADR-0011](0011-history-and-naming.md), [ADR-0020](0020-cpp-c-abi-and-python-bindings.md), [ADR-0021](0021-shape-derivatives.md), [ADR-0026](0026-modelling-operations-by-contracts.md), [ADR-0027](0027-controlled-topology-mutation-and-euler-operators.md), [ADR-0031](0031-attribute-propagation-and-conflicts.md).

This proposal supplements the existing history decision and supersedes none. All listed verification is required as the corresponding operations become supported; no history implementation is claimed here.
