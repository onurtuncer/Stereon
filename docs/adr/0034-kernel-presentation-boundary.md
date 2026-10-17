# ADR-0034: Separate the Modelling Kernel from Presentation and Interaction

- **Status:** Proposed
- **Date:** 2026-10-16
- **Owner:** Onur Tuncer
- **Phase:** 0 (boundary); applied as viewing, selection and annotation facilities are introduced
- **Scope:** The boundary between the modelling kernel and presentation, interaction, text and display caches
- **Related:** ADR-0004, ADR-0005, ADR-0007, ADR-0014, ADR-0028, ADR-0030, ADR-0033

## Context

Stereon will need interactive inspection of curves, surfaces, topology, and modelling results, followed by selection, dimensions, annotations, and drawing workflows. These facilities must not make the modelling kernel dependent on a graphics backend, window system, font engine, or interactive session.

A displayed object is a representation of model data. It is not the model entity itself. A single entity may have several representations in different views, and a representation may aggregate several entities. Display tessellation and picking results are approximations whose validity depends on the model revision and presentation settings.

This decision defines the boundary now without requiring a complete interactive framework in the initial kernel implementation.

## Options

1. **Put rendering and selection state on kernel entities.** Couples model lifetime and algorithms to frontend resources and complicates headless use and multiple views.
2. **Use rendered primitives as model identity.** Remeshing, backend changes, and topology edits invalidate those identifiers.
3. **Build a complete interaction framework immediately.** The boundary can be enforced before advanced viewer features exist.
4. **Use external viewers indefinitely.** Useful for early inspection, but insufficient for integrated entity selection and modelling diagnostics.
5. **An explicit boundary, with presentation and interaction as adapters over public kernel interfaces, delivered incrementally.**

## Decision

Adopt option 5. Options 1 and 2 are rejected; option 3 is deferred; option 4 serves only early inspection.

### 1. Dependencies point toward the kernel

The modelling kernel owns geometry, topology, modelling algorithms, numerical contracts, and operation results. It MUST build and run headlessly without a renderer, windowing toolkit, font library, GPU context, or event loop.

Presentation adapters MAY depend on public kernel interfaces. The kernel MUST NOT import or link presentation, interaction, GUI, or rendering modules. Public kernel interfaces MUST NOT expose backend-specific types, graphics handles, presentation objects, or UI callbacks.

An application/controller layer coordinates modelling commands and presentation updates. User input reaches modelling operations through this layer; rendering and selection code MUST NOT mutate kernel storage directly.

### 2. Ownership is explicit

| Concern | Owner |
| --- | --- |
| Geometry, topology, model validity, operation outcomes | Modelling kernel |
| Authoritative entity identity and reference validation | Model/kernel API |
| Saved annotations, units, dimension semantics, model-linked styles | Document/application layer |
| Display objects, cameras, visibility overrides, highlighting, diagnostic labels | Presentation/interaction layer |
| Display tessellation caches, glyph caches, GPU buffers, shaders | Presentation/backend layer |
| Font loading, text layout, glyph rendering | Optional text service used by presentation/export |

Persistent material or colour attributes may belong to a model/document. Their visual interpretation and temporary overrides belong to presentation. This ADR does not prohibit reusable headless meshing algorithms in the kernel: presentation owns the display-specific mesh cache, not necessarily the algorithm that produces it.

### 3. Selection resolves to validated model references

Picking MUST return a backend-independent selection result containing the model/document identity, entity reference, and model revision or equivalent validity token. Assembly instance identity and subentity identity MUST be included where needed to distinguish occurrences.

GPU object IDs, triangle indices, raw addresses, and unqualified storage indices MUST NOT serve as authoritative model references. A rendering adapter maintains the mapping from its primitives to model references. The model validates references before an operation uses them; stale or deleted references MUST produce an explicit failure or require a fresh selection.

Presentation identifiers are independent of model identifiers. Model references do not imply persistent naming across edits. Splits, merges, deletions, and replacements must be handled through operation results and lineage where available; ambiguous mappings MUST NOT silently select an arbitrary descendant.

### 4. Presentations are derived and revision-aware

Display representations MUST be rebuildable from model/document data and presentation settings. Cache validity MUST account for the relevant model revision, instance placement, and rendering or tessellation parameters.

The application passes committed changes, operation results, or revision notifications to presentation. Kernel operations MUST NOT require a viewer to consume notifications or acknowledge updates. A full presentation rebuild is an acceptable initial implementation; incremental invalidation is a later optimisation.

Background presentation work MUST use an immutable snapshot or another explicitly documented safe read mechanism. A result computed for an obsolete revision must be discarded or deliberately retained as a visibly separate snapshot; it must not silently replace the current representation.

### 5. Display accuracy is distinct from modelling accuracy

Picking and display tessellation tolerances MUST remain separate from modelling tolerances. Changing zoom, display quality, or GPU precision MUST NOT alter model geometry or modelling outcomes.

A triangle hit can identify a candidate entity and approximate hit point. Commands requiring an exact geometric point MUST refine or validate that point through kernel geometry queries with explicit numerical requirements.

### 6. Text and dimensions retain their meaning

Ordinary labels and annotations MUST remain text plus style and placement data. They MUST NOT require conversion into B-rep geometry for display. Font and graphics resources remain outside the kernel.

Associative dimensions MUST retain their geometric references, measurement definition, units, and formatting separately from their displayed string. Invalidated references must make the dimension explicitly unresolved rather than leave a stale value appearing authoritative. Purely presentational notes are permitted and must not masquerade as associative measurements.

Text-to-geometry conversion is an explicit optional operation in an adapter/service. It may pass generated curves or profiles to kernel modelling APIs for engraving or embossing. The kernel operates on that geometry without needing to know the source text or font engine.

### 7. Delivery is incremental

1. During kernel foundations, use headless tests and export-based inspection.
2. Alongside the first meaningful B-rep operations, introduce a thin inspection viewer with navigation, edge/face display, and diagnostic labels.
3. Add validated picking, highlighting, selection filters, and revision-driven updates when interactive modelling requires them.
4. Add associative dimensions, drawing presentations, and richer annotation behaviour when their document workflows are implemented.

This ADR does not choose a graphics API, UI toolkit, font library, module naming scheme, or drawing file format. Those choices may be recorded independently. It does not require cloning OCCT's AIS API or inheritance structure.

## Consequences

- The kernel remains usable in tests, servers, batch processing, bindings, and alternative frontends. Rendering backends can evolve without changing modelling contracts, and several views can share one model safely.
- The cost is an explicit adapter layer, model-to-presentation mappings, and revision/invalidation handling. Interactive prototypes must respect those boundaries even if a direct shortcut initially looks simpler.

## Verification

- CI can build and exercise the kernel with all presentation and GUI components disabled.
- Build dependency checks reject reverse dependencies from the kernel to presentation/backend modules.
- When selection is introduced, stale references and distinct assembly occurrences are covered by meaningful integration checks.
- Changing display tessellation does not change model identity or modelling results.
- When asynchronous presentation is introduced, obsolete results cannot overwrite the current revision's presentation.
- Ordinary annotation display does not create B-rep entities; text-to-geometry is an explicit action.

These criteria apply as the corresponding facilities are introduced; they do not require implementing later-stage features now.

## Relationship to existing decisions

- [ADR-0033](0033-modular-sdk-product-scope.md), the companion decision, defines the broader SDK and product scope. This ADR details the kernel–presentation boundary within it.
- [ADR-0004](0004-topology-storage-and-handles.md): typed handles carry a store and a generation, and stale handles are detected. These are the kernel-side validity tokens that selection results resolve to (section 3). Handles remain local to one store, so cross-edit identity comes from lineage, not handles.
- [ADR-0005](0005-immutability-and-sharing.md) and [ADR-0007](0007-threading-and-determinism.md): immutable shapes and snapshots are the safe read mechanism for background presentation work (section 4), and presentation never alters a deterministic modelling result.
- [ADR-0014](0014-mesh-generation-boundary.md) keeps display tessellation algorithms in the kernel's family. This ADR keeps the display-specific mesh cache in presentation (section 2), consistent with it.
- [ADR-0028](0028-purpose-specific-numerical-tolerances.md): display and picking tolerances are a purpose of their own, never interchangeable with modelling tolerances (section 5).
- [ADR-0030](0030-entity-lineage-in-operation-results.md): operation lineage is how selections and associative dimensions follow splits, merges and replacements, without silently picking an arbitrary descendant (sections 3 and 6).
- In the current repository, no library under `src/` depends on a graphics, window or font library, and `tools/check_layering.py` already rejects upward dependencies between libraries. The *build dependency check* criterion extends it to presentation and backend modules when they are added. The Phase 0 debug viewer (TODO.md §8) is the first consumer bound by this boundary.

This proposal supplements these decisions and supersedes none.
