# ADR-0011: Operation history and topological naming

- **Status:** Proposed (to be accepted before Phase 3)
- **Date:** 2026-09-28
- **Revised:** 2026-10-10 (`OperationResult`; remap from `Compact()`; attribute conflict policy)
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

- Every shape-producing operation returns a `History` inside its `OperationResult` (ADR-0006) that maps each input entity to output entities with a tag: `Generated`, `Modified`, `Deleted` or `Unchanged`.
- `Freeze()` keeps indices (ADR-0004), so no remap is needed for it. The remap table produced by `Compact()` is composed into the history automatically, so algorithms only record their own semantic relations.
- Generated entities record their generators (e.g. a fillet face records the edge it replaced and the two faces it joins; the wall of a drilled hole records the tool's cylindrical face).
- **Attribute propagation policy.** Sparse attributes (ADR-0004) are propagated through history by default: `Modified` and `Unchanged` entities keep their attributes, and `Generated` entities inherit from their generators of the same entity kind. When an output entity has several ancestors with different values for one attribute key (two tagged faces merged into one, a hole wall generated from a tagged tool face into a tagged body), the kernel applies the `MergePolicy` for that key from the `Context`: `KeepFirst` (lowest input index, deterministic), `Drop`, or a user callback that receives the candidate values and the entities. Every conflict, however resolved, is reported as a warning in the `OperationResult`, so an application can refuse a result whose boundary-condition tags were decided by default.
- A persistent-ID attribute slot is reserved but not interpreted by the kernel in v1.0.

## Consequences

- Applications get reliable attribute propagation without geometric guessing.
- Every algorithm must record history; review and tests enforce it.
- A future naming system can be built from history plus generator records without changing the core.

## Verification

For every constructive operation and boolean in the corpus, check that each output entity is reachable from at least one input or is tagged `Generated` with generators, and that attributes survive `Modified` and `Unchanged` relations. The duct example (ADR-0014): a duct body with faces tagged `inlet`, `outlet` and `wall`, cut by a cylinder whose faces are tagged `wall`, yields a result in which every face carries exactly one tag, the hole wall is `wall`, the inlet and outlet keep their tags, and the `OperationResult` lists every conflict that was resolved by policy.
