# ADR-0004: Topology storage and handles

- **Status:** Proposed
- **Date:** 2026-09-28
- **Owner:** Onur Tuncer
- **Phase:** 0

## Context

Booleans constantly ask adjacency questions ("which faces share this edge?", "which edges meet at this vertex?"). The topology store must make these fast and cache-friendly, support the immutable snapshots of ADR-0005, give entities stable identity for history maps (ADR-0011) and the native format (ADR-0012), and represent non-manifold intermediate results (sheets, wires, edges shared by three or more faces).

## Options

1. **Pointer graph of reference-counted objects** (OCCT `TopoDS_TShape`). Familiar; poor cache locality, atomic reference-count contention, identity by memory address, hard to serialise.
2. **Structure-of-arrays per entity type with typed index handles.** Cache-friendly, trivially serialisable, cheap to copy, thread-safe when frozen; raw indices are opaque in a debugger.
3. **ECS-style (entities plus optional component tables).** Arbitrary attributes; topology traversal becomes a chain of table lookups.

## Decision

Adopt option 2, with a small ECS-style side table for sparse attributes.

1. **Builder and frozen shape.** `ShapeBuilder` is mutable and single-threaded; its handles are 32-bit index + 32-bit generation with free lists, so stale handles are detected. `freeze()` compacts to an immutable `Shape` with plain 32-bit indices and emits an old→new remap table that feeds the history maps.
2. **Coedge core with radial rings.** Each coedge stores `next`, `prev`, `radial` (next use of the same edge, forming a ring), `edge`, `loop`, `pcurve` and `reversed`. The radial ring handles manifold solids, sheets and non-manifold edges with one structure.
3. **Geometry in separate pools.** Topology references curves, pcurves and surfaces by index into pools of type-erased values (ADR-0003); several faces may share one surface.
4. **No per-sub-shape transforms.** Bodies are stored in their own coordinates. Instancing belongs to an assembly layer above the kernel.
5. **Copy-on-write pages.** Arrays are stored in 256-entry pages behind `shared_ptr<const Page>`; an operation copies only the pages it touches (ADR-0005).
6. **Sparse attributes** (names, colours, user tags, persistent IDs) live in a side table keyed by handle, keeping hot arrays small.

Handles are local to one `Shape`. Identity across operations comes from history maps, never from handles.

## Consequences

- Parallel read-only algorithms, serialisation, determinism and shape diffing in tests become straightforward.
- Debugging needs pretty-printers and the debug viewer from Phase 0.
- We must build adjacency range views (`faces_of(edge)`, `edges_of(vertex)`, loop walks), the remap machinery, and an array-level validity checker.
- Non-manifold support from day 1 costs some complexity in Euler operators but avoids rewriting the boolean code later.

## Verification

Phase 0 prototype with a ~10⁶-face body. Pass: adjacency queries ≥ 5× faster than OCCT `TopExp::MapShapesAndAncestors`; memory per face < 50% of OCCT; modifying one face copies < 1% of pages.
