# ADR-0011: Operation history and topological naming

- **Status:** Proposed (to be accepted before Phase 3)
- **Date:** 2026-09-28
- **Revised:** 2026-10-10 (`OperationResult`; remap from `Compact()`; attribute conflict policy)
- **Revised:** 2026-10-14 (qualified lineage and per-key propagation extensions in ADR-0030/0031)
- **Owner:** Onur Tuncer
- **Phase:** 3

## Context

Applications built on a kernel need to know what happened to each entity during an operation: which face of the result came from which input face, which edges are new. This is needed for attributes (colours, names), for downstream references (a fillet applied to "this edge") and, eventually, for parametric modelling. The topological naming problem — keeping references stable when a model is rebuilt with new parameters — is notoriously hard, and parametric history is out of scope for v1.0.

## Options

1. **No history.** Simplest; pushes all tracking onto applications, which then guess by geometry.
2. **History maps per operation.** Each operation reports input→output relations; applications build their own naming on top.
3. **Full persistent naming inside the kernel.** Most capable; large scope and a research problem in its own right.

## Decision

Adopt option 2 now, designed so option 3 can be built on it later.

- Every shape-producing operation returns a `History` inside its `OperationResult` (ADR-0006). [ADR-0030](0030-entity-lineage-in-operation-results.md) defines qualified input/output references, many-to-many continuation and generation relations, Deleted dispositions, explicit split/merge event groups, creation roots for primitives/import, composition and atomic publication. The relation vocabulary remains `Generated`, `Modified`, `Deleted` and `Unchanged`; Deleted has no output target and may coexist with generation provenance.
- `Freeze()` keeps indices (ADR-0004), so no remap is needed for it. The remap table produced by `Compact()` is composed into the history automatically, so algorithms only record their own semantic relations.
- Generated entities record their generators (e.g. a fillet face records the edge it replaced and the two faces it joins; the wall of a drilled hole records the tool's cylindrical face). Operations without input topology record an operation/import creation root, not a fictitious source entity. Generation is distinct from survival or a split.
- **Attribute propagation policy.** [ADR-0031](0031-attribute-propagation-and-conflicts.md) defines the detailed per-key transfer contract. Passive same-kind continuations and Generated outputs inherit eligible attributes; semantic keys use registered split/merge, revalidation and required-coverage policies. For conflicting candidates, `MergePolicy` offers `KeepFirst` (ordered by input ordinal, kind, index and stable role/event key), `Drop`, or a deterministic custom resolver. Conflicts and nontrivial losses are structured warnings in `OperationResult`; strict policies promote them to failure before publication. Geometry-dependent data, conserved quantities and unique IDs are not generically guaranteed correct by passive copying.
- A persistent-ID attribute slot is reserved but not interpreted by the kernel in v1.0.

## Consequences

- Applications get reliable attribute propagation without geometric guessing.
- Every algorithm must record history; review and tests enforce it.
- A future naming system can be built from history plus generator records without changing the core.

## Verification

For every constructive operation and boolean in the corpus, check that each output entity is covered by a qualified continuation or Generated provenance (input generators or an operation/import creation root), and that passive attributes survive eligible `Modified` and `Unchanged` relations, while registered semantic keys satisfy ADR-0031's transfer/revalidation rules. Exercise ADR-0030 split/merge composition, deletion and empty-result histories as well. The duct example (ADR-0014): a duct body with faces tagged `inlet`, `outlet` and `wall`, cut by a cylinder whose faces are tagged `wall`, yields a result in which every face carries exactly one tag, the hole wall is `wall`, the inlet and outlet keep their tags, and the `OperationResult` lists every conflict that was resolved by policy.
