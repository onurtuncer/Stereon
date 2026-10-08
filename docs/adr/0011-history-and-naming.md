# ADR-0011: Operation history and topological naming

- **Status:** Proposed (to be accepted before Phase 3)
- **Date:** 2026-09-28
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

- Every operation returns, alongside the result, a `History` that maps each input entity to output entities with a tag: `Generated`, `Modified`, `Deleted` or `Unchanged`.
- The remap table produced by `ShapeBuilder::Freeze()` (ADR-0004) is composed into the history automatically, so algorithms only record their own semantic relations.
- Generated entities record their generators (e.g. a fillet face records the edge it replaced and the two faces it joins).
- Sparse attributes (ADR-0004) are propagated through history by default: `Modified` and `Unchanged` entities keep their attributes.
- A persistent-ID attribute slot is reserved but not interpreted by the kernel in v1.0.

## Consequences

- Applications get reliable attribute propagation without geometric guessing.
- Every algorithm must record history; review and tests enforce it.
- A future naming system can be built from history plus generator records without changing the core.

## Verification

For every constructive operation and boolean in the corpus, check that each output entity is reachable from at least one input or is tagged `Generated` with generators, and that attributes survive `Modified` and `Unchanged` relations.
