# ADR-0004: Topology storage and handles

- **Status:** Proposed
- **Date:** 2026-10-08
- **Owner:** Onur Tuncer
- **Phase:** 0

## Context

Booleans constantly ask adjacency questions ("which faces share this edge?", "which edges meet at this vertex?"). The topology store must make these fast and cache-friendly, support the immutable snapshots of ADR-0005, give entities stable identity for history maps (ADR-0011) and the native format (ADR-0012), and represent non-manifold intermediate results (sheets, wires, edges shared by three or more faces).

This ADR separates the **logical model** (which entities exist and how they relate) from the **physical storage** (how they sit in memory). The logical model follows Golovanov (*Geometric Modeling: The Mathematics of Shapes*, the design basis of the C3D kernel): solid → shell → face → loop → oriented edge (coedge) → edge → vertex, with explicit orientation at each level and edges defined by the surfaces they join (ADR-0009). The physical storage below is Stereon's own; C3D's pointer-linked, reference-counted objects are not used, and no C3D code or API is copied.

### Prior art: OpenCASCADE `BRepGraph`

As of OCCT master (August 2026), OpenCASCADE contains `BRepGraph` (`src/ModelingData/TKBRep/BRepGraph`, backed by `BRepGraphInc`), an index-based incidence-table representation of B-rep topology. It independently arrives at much of this ADR:

- per-kind entity tables with typed 32-bit IDs (`BRepGraph_NodeId::Typed<Kind>`, cross-kind conversion deleted);
- a `CoEdge` entity owning the pcurve of each edge–face use (half-edge pattern);
- generation counters, soft removal, compaction, per-kind UIDs, a history layer, a validation pass, fuzz tests and a parallel policy;
- assembly tables (`Product`, `Occurrence`) separate from part topology.

It differs from Stereon in ways that matter for this decision:

| | OCCT `BRepGraph` | Stereon |
| --- | --- | --- |
| Role | An additional layer, populated from and reconstructed to `TopoDS`; the boolean and fillet toolkits (`TKBO`, `TKFillet`) do not use it in the inspected snapshot | The only topology store; every algorithm works on it natively |
| Mutability | A mutable graph with an editor, version stamps and cache invalidation | Mutable `ShapeBuilder`, then an immutable frozen `Shape` with copy-on-write pages (ADR-0005) |
| Geometry | Representation tables holding mutable, reference-counted `Geom_*` handles | Immutable type-erased geometry values in pools (ADR-0003) |
| Edge curves | Independent 3D curve plus pcurves, consistency via SameParameter/SameRange semantics | Edges defined by their two surfaces and pcurves on a shared parameter; 3D curve is a cache (ADR-0009) |
| Placement | `LocalLocation` on child and occurrence references | No transforms inside a body; placement only in the assembly layer |
| Non-manifold | Edge→coedge relation lists | Radial coedge rings |

**Implications.** `BRepGraph` confirms the direction of this ADR, so the storage layout alone is not a differentiator; Stereon's case rests on native use throughout, immutability with determinism, the edge model and the tolerance discipline. Its relation tables, history layer, deduplication, validation and fuzz tests are worth studying. Its tables also offer a convenient bridge for the OCCT differential-testing harness. OCCT is LGPL-2.1: its design is studied for ideas only, and no OCCT code is copied into Stereon (ADR-0008). If OCCT moves its core algorithms onto `BRepGraph`, revisit Stereon's positioning in the roadmap.

## Options

1. **Pointer graph of reference-counted objects** (OCCT `TopoDS_TShape`, C3D). Familiar; poor cache locality, atomic reference-count contention, identity by memory address, hard to serialise.
2. **Structure-of-arrays per entity type with typed index handles.** Cache-friendly, trivially serialisable, cheap to copy, thread-safe when frozen; raw indices are opaque in a debugger.
3. **ECS-style (entities plus optional component tables).** Arbitrary attributes; topology traversal becomes a chain of table lookups.

## Decision

Adopt option 2, with a small ECS-style side table for sparse attributes.

1. **Builder and frozen shape.** `ShapeBuilder` is mutable and single-threaded; its handles are 32-bit index + 32-bit generation with free lists, so stale handles are detected. `Freeze()` compacts to an immutable `Shape` with plain 32-bit indices and emits an old→new remap table that feeds the history maps.
2. **Coedge core with radial rings.** Each coedge stores `Next`, `Prev`, `Radial` (next use of the same edge, forming a ring), `Edge`, `Loop`, `PCurve` and `Reversed`. The radial ring handles manifold solids, sheets and non-manifold edges with one structure.
3. **Geometry in separate pools.** Topology references curves, pcurves and surfaces by index into pools of type-erased values (ADR-0003); several faces may share one surface.
   **Edge records follow ADR-0009.** Each edge stores its curve kind (`Intersection`, `Exact`, `Boundary`, `Degenerate`), its tolerance, and either its exact 3D curve (`Exact`) or a cached 3D approximation with its certified deviation (`Intersection`). The pcurves stay on the coedges, all on the edge's shared parameter; the surfaces are reached through each coedge's face, so an edge never stores its own surface references. Cached 3D curves are invalidated whenever a pcurve of the edge changes.
4. **No per-sub-shape transforms.** Bodies are stored in their own coordinates. Instancing belongs to an assembly layer above the kernel.
5. **Copy-on-write pages.** Arrays are stored in 256-entry pages behind `shared_ptr<const Page>`; an operation copies only the pages it touches (ADR-0005).
6. **Sparse attributes** (names, colours, user tags, persistent IDs) live in a side table keyed by handle, keeping hot arrays small.

Handles are local to one `Shape`. Identity across operations comes from history maps, never from handles.

### Debuggability safeguards

Index handles fail differently from pointers: a stale or off-by-one index usually points to a *valid but wrong* entity and produces a plausible, incorrect result instead of a crash. The following safeguards are part of this decision and must exist before Phase 1 algorithms are written.

1. **Typed handles.** `VertexId`, `EdgeId`, `CoedgeId`, `LoopId`, `FaceId`, `ShellId`, `SolidId` are distinct types (`Id<Tag>`). Mixing entity kinds is a compile error; there are no implicit conversions to or from integers.
2. **Checked handles in debug builds.** With `STEREON_CHECKED_HANDLES` (on in `debug`, `asan` and fuzzing presets), every handle also carries the 32-bit ID of the store it came from and its generation. Every access verifies store ID, bounds and generation, and a mismatch is a contract violation reporting the handle, the expected store and the entity kind. Release builds compile handles down to a plain 32-bit index; `sizeof` checks in CI confirm this.
3. **Orientation bits are explicit.** Where an orientation flag is packed into a handle (e.g. the sense bit of an oriented coedge), it is accessed only through named functions (`IsReversed(h)`, `BaseOf(h)`), never by manual masking.
4. **Pretty-printers show the entity, not the index.** GDB (Python) and LLDB formatters resolve a handle against its store and display it, e.g. `Edge#4711 [Line, tol 1e-5, faces 12|37]`. A `Describe(shape, handle)` function returns the same text for logs and assertions.
5. **Validity checking after every operation.** In debug builds, every public operation runs the array-level validity checker on its result, so corruption is reported at the operation that caused it. The checker names the offending entities with `Describe()`.
6. **Replayable journals.** Each operation can record its inputs and parameters to a journal file; any failure can be replayed deterministically in isolation (ADR-0007).
7. **Visual and textual inspection.** Any shape or sub-shape can be sent to the debug viewer (`STEREON_DUMP(obj)`) and written as JSON with `stn-dump` (ADR-0012), so two versions of a shape can be diffed.

**Fallback.** If the Phase 0 prototype shows debugging remains too costly despite these safeguards, switch to a hybrid: pointer-linked entities in the mutable `ShapeBuilder`, indices only in the frozen `Shape`. This must be decided before Phase 1 algorithms start, when the change costs weeks rather than a rewrite.

## Consequences

- Parallel read-only algorithms, serialisation, determinism and shape diffing in tests become straightforward.
- Debugging depends on the safeguards above; they are Phase 0 deliverables, not later polish.
- Debug builds are slower and handles are larger (checked handles plus per-operation validation); acceptable because release builds pay nothing.
- We must build adjacency range views (`FacesOf(edge)`, `EdgesOf(vertex)`, loop walks), the remap machinery, and an array-level validity checker.
- Non-manifold support from day 1 costs some complexity in Euler operators but avoids rewriting the boolean code later.

## Verification

Phase 0 prototype with a ~10⁶-face body. Pass:

- Adjacency queries ≥ 5× faster than OCCT `TopExp::MapShapesAndAncestors`; memory per face < 50% of OCCT; modifying one face copies < 1% of pages.
- Release-build handles are exactly 4 bytes (`static_assert` + CI check).
- A seeded-fault test suite — stale handles, handles from another shape, off-by-one indices, wrong entity kinds, corrupted adjacency — is caught in debug builds 100% of the time, each at the first faulty access or at the end of the operation that caused it, with a `Describe()` message naming the entity.
- Pretty-printers render every handle type in GDB and LLDB (tested in CI with scripted debugger sessions).
