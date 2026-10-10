# ADR-0005: Immutability and structural sharing

- **Status:** Proposed
- **Date:** 2026-09-28
- **Revised:** 2026-10-10 (sharing requires index-stable `Freeze()`; `OperationResult`; page size measured, not fixed)
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

- `Shape` is immutable after `Freeze()` (ADR-0004). All kernel operations take `const Shape&` and return `std::expected<OperationResult, KernelError>`, where `OperationResult` carries the new `Shape` together with its history and diagnostics (ADR-0006).
- Storage is paged; pages are `shared_ptr<const Page>`. Builders copy a page on first write. Sharing only works because `Freeze()` keeps indices stable and leaves tombstones for deleted entities (ADR-0004); renumbering is the explicit `Compact()` step, which by design copies everything.
- The page size per table is a tuning constant, initially 256 entries. The verification below measures bytes per page, share ratio and copy cost on representative workloads for 128, 256 and 1 024 entries and fixes the values; nothing in the public API depends on them.
- Geometry pools follow the same scheme; large NURBS payloads are shared immutable buffers.
- Undo is keeping the previous `Shape` value. No kernel-level undo log exists.

## Consequences

- Shapes can be shared across threads without locks; ADR-0007 relies on this.
- Memory use stays proportional to what changed, not to model size.
- Reference-count traffic on pages must stay off hot paths: algorithms take raw spans of a page once, not per entity.
- In-place editing APIs do not exist; users who want them build on `ShapeBuilder`.

## Verification

Benchmark: a sequence of 100 local edits on a 10⁶-face body, mixing surface replacement, face deletion, face splitting and edge splitting in equal parts, with no `Compact()` in the sequence. Pass: total memory < 1.2× the single body; each edit's copy cost < 1% of the body's storage. A `Compact()` at the end returns the body to within 1.05× of a freshly built one.
