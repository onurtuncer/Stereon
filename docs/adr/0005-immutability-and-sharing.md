# ADR-0005: Immutability and structural sharing

- **Status:** Proposed
- **Date:** 2026-09-28
- **Owner:** Onur Tuncer
- **Phase:** 0

## Context

CAD applications need undo, cheap copies of large bodies, and safe concurrent reads (tessellation on one thread while the UI queries another). Mutable shapes with locks or undo logs make all three error-prone.

## Options

1. **Mutable shapes with an undo log.** Simple for single-threaded use; concurrency needs locking, and undo logs grow complex with every new operation.
2. **Deep copy on every operation.** Simple and safe; far too expensive for large bodies.
3. **Immutable persistent shapes with copy-on-write pages.** Operations return a new `Shape` that shares unchanged storage with its input.

## Decision

Adopt option 3.

- `Shape` is immutable after `freeze()` (ADR-0004). All kernel operations take `const Shape&` and return `std::expected<Shape, KernelError>`.
- Storage is paged (256 entries per page); pages are `shared_ptr<const Page>`. Builders copy a page on first write.
- Geometry pools follow the same scheme; large NURBS payloads are shared immutable buffers.
- Undo is keeping the previous `Shape` value. No kernel-level undo log exists.

## Consequences

- Shapes can be shared across threads without locks; ADR-0007 relies on this.
- Memory use stays proportional to what changed, not to model size.
- Reference-count traffic on pages must stay off hot paths: algorithms take raw spans of a page once, not per entity.
- In-place editing APIs do not exist; users who want them build on `ShapeBuilder`.

## Verification

Benchmark: a sequence of 100 local edits on a 10⁶-face body. Pass: total memory < 1.2× the single body; each edit's copy cost < 1% of the body's storage.
