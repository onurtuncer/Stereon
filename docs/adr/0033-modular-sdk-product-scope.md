# ADR-0033: Structure Stereon as a Modular CAD SDK with an Independent Kernel

- **Status:** Proposed
- **Date:** 2026-10-16
- **Owner:** Onur Tuncer
- **Phase:** 0 (scope and dependency policy); applied as each capability is introduced
- **Scope:** The product scope of Stereon and the responsibility and dependency boundaries between its capabilities
- **Related:** ADR-0012, ADR-0014, ADR-0019, ADR-0020, ADR-0030, ADR-0031, ADR-0034

## Context

Stereon's ambition extends beyond geometric algorithms to a reusable CAD development platform. Applications need model organisation, meshing, exchange, visualisation, interaction, and eventually annotations and drawings. Treating all these capabilities as one indivisible kernel would introduce unnecessary dependencies and make headless deployment difficult. Conversely, leaving them entirely to unrelated applications would duplicate CAD-specific integration and weaken interoperability within Stereon.

We need an explicit product scope and a division of responsibilities before these capabilities grow. This ADR defines target architectural responsibilities, not an inventory of existing repository modules or a commitment to implement every capability immediately.

## Options

1. **Deliver only a kernel indefinitely.** Narrower scope, but leaves recurring CAD integration work to each application.
2. **Build one monolithic CAD application.** Provides a direct user experience but couples reusable facilities to product workflows and frontend dependencies.
3. **Develop the whole SDK before validating the kernel.** Spreads effort across speculative abstractions and delays useful feedback.
4. **Adopt an existing platform's class hierarchy wholesale.** May ease familiarity but would constrain Stereon without proving that the design suits its contracts.
5. **A modular SDK around an independently usable kernel, delivered incrementally.**

## Decision

Adopt option 5. Option 1 is rejected as the long-term platform scope but retained as a supported standalone deliverable (the kernel). Options 2–4 are rejected; no compatibility commitment to another platform's hierarchy is made.

### 1. Distinguish three deliverables

Stereon SHALL be developed as a modular CAD SDK with three distinct deliverables:

1. **Stereon kernel:** independently usable geometry, topology, numerical foundations, and modelling operations. It supports headless execution.
2. **Stereon SDK:** the kernel together with independently selectable supporting capabilities, subject to their declared dependencies.
3. **Stereon inspection application:** a small application composed from SDK modules for inspecting results, diagnosing defects, and demonstrating workflows.

A full end-user CAD editor is outside the scope of this decision. Parametric feature authoring, sketch-solving workflows, manufacturing tools, and a complete drawing editor require separate scope decisions. The inspection application MUST NOT become an implicit commitment to those products.

### 2. Establish responsibility boundaries

The following are logical capability boundaries. They do not prescribe exact library names, C++ module names, package counts, or process boundaries.

| Capability | Responsibility | Dependency rule |
| --- | --- | --- |
| Foundations and modelling kernel | Mathematical primitives, numerical contracts, geometry, topology, intersections, modelling operations, and operation results | No dependencies on documents, presentation, GUI, or graphics backends |
| Meshing | Reusable headless tessellation and mesh-generation facilities | May depend on kernel APIs; must not require a viewer |
| Model/document services | Parts, assembly instances, attributes, document references, and persistent document state | May depend on the kernel; must remain headless |
| Data exchange | Import/export, format mappings, and reports of unsupported or degraded information | Geometry adapters may depend on the kernel; richer adapters may also depend on document services |
| Visualisation | Presentations, cameras, display styles, diagnostic labels, and rendering integration | Depends on kernel/meshing interfaces as needed; document integration is an optional adapter |
| Interaction | Picking, highlighting, selection filters, and manipulators | Depends on presentation and validated model-reference contracts; delegates edits to application services |
| Application services | Command coordination, transactions, undo/redo, and propagation of committed changes | Coordinates headless modelling/document services; UI-specific controllers sit above it |
| Annotations and drawings | Associative measurement definitions, annotation data, projected-view generation, sheet data, and their presentations | Headless semantics and computation are separated from optional interactive presentation |
| Language/API bindings | Expose selected public capabilities to host applications and languages | Depend on exposed modules; internal modules must not depend on bindings |

Core transactional guarantees belong to the service that owns the affected state. Application services coordinate those guarantees across commands; placing undo/redo here does not permit partially applied or invalid kernel operations.

### 3. Keep the dependency graph acyclic and capabilities optional

Each module MUST declare its dependencies and public contracts. Higher-level integration belongs in adapters or composition layers rather than reverse dependencies in lower-level modules.

The kernel MUST build and run without the inspection application, presentation, interaction, font services, windowing libraries, or a GPU. Headless document, exchange, meshing, and annotation workflows MUST likewise avoid mandatory UI dependencies.

The SDK is a family of compatible capabilities, not a requirement to link every module. Enabling a capability may enable its declared dependencies; it must not silently pull in unrelated application or GUI features.

Shared public data types should live in the lowest appropriate owning module. A general-purpose shared module MUST NOT become a dumping ground used to disguise dependency cycles.

### 4. Own CAD contracts and integrate supporting technologies

Stereon owns its public CAD contracts, entity-reference semantics, operation outcomes, and adaptation between modules. It need not implement every supporting technology itself.

Existing libraries may provide rendering, font processing, window management, UI widgets, and other infrastructure. Such choices require evaluation of suitability, licensing (ADR-0008), deployment, and maintenance. This ADR selects no particular dependency or backend.

Third-party implementation types SHOULD remain behind adapters where exposing them would unnecessarily constrain Stereon's public interfaces. Public exposure, where justified, must be an explicit architectural choice.

### 5. Keep the inspection application thin

The inspection application SHALL consume public SDK interfaces. Reusable modelling, meshing, exchange, and document behaviour MUST live in the appropriate SDK module rather than inside the application.

Its initial remit is loading or constructing test models, navigating views, inspecting entities, showing diagnostics, and exercising selected operations. It is a consumer and integration check of the SDK, not the owner of modelling semantics.

Interactive inspection complements automated validation. Visual appearance alone does not establish geometric or topological correctness.

### 6. Deliver capabilities incrementally

The immediate priority is a sound kernel and a thin inspection path. Supporting modules are introduced when a concrete workflow requires them:

- Establish kernel contracts and headless validation first.
- Introduce reusable tessellation and a thin viewer alongside meaningful B-rep inspection needs; export-based inspection may serve earlier work.
- Add document services, validated interaction, and command coordination as editing and assembly workflows require them.
- Expand exchange and bindings incrementally around stable contracts and specific interoperability needs.
- Introduce annotation and drawing capabilities after the necessary geometry, reference, and document contracts exist.

This is a dependency-informed delivery policy, not a dated roadmap or a replacement for existing milestones. It neither requires empty placeholder packages nor delays all exchange or binding work until the kernel is complete.

## Consequences

- Applications can adopt the kernel alone or select a larger set of compatible CAD capabilities. Headless services and alternative frontends remain viable. The inspection application can evolve without defining the kernel's runtime environment.
- The SDK ambition increases the long-term workload beyond kernel development. Explicit module contracts, integration testing, dependency management, and documentation are required. Architectural separation does not itself deliver feature parity with another CAD platform.
- Module boundaries may be refined as implementations reveal actual coupling. Moving responsibilities across these boundaries must preserve the dependency rules or be justified by a superseding ADR.

## Verification

- The build supports a kernel-only configuration without presentation or GUI dependencies.
- Modules have documented owners, public interfaces, and an acyclic declared dependency graph.
- New modules are justified by concrete consumers or workflows.
- The inspection application uses public SDK interfaces; reusable algorithms do not reside only in its source tree.
- Headless services can be exercised without creating a window or rendering context.
- Integration checks are added at actual module boundaries as those modules appear, rather than creating tests for unimplemented capabilities.

## Relationship to existing decisions

- [ADR-0034](0034-kernel-presentation-boundary.md), the companion decision, specifies presentation ownership, model selection references, revision handling, and text/display behaviour in greater detail. This ADR defines the broader SDK and product scope and does not supersede it.
- [ADR-0014](0014-mesh-generation-boundary.md) places display and analysis-grade surface meshing in the kernel's family and volume meshing in adapters. That is the *Meshing* capability above: headless, kernel-dependent, never requiring a viewer.
- [ADR-0019](0019-library-linking-packaging-and-abi.md) and [ADR-0020](0020-cpp-c-abi-and-python-bindings.md) govern how the deliverables are packaged and exposed. Their aggregate target and staged C ABI are where *independently selectable* capabilities and *bindings* become concrete.
- [ADR-0012](0012-native-file-format.md) (native persistence) is the first *Data exchange* capability. Document state beyond geometry and shapes belongs to *Model/document services*.
- [ADR-0030](0030-entity-lineage-in-operation-results.md) and [ADR-0031](0031-attribute-propagation-and-conflicts.md) define the lineage and attribute semantics that document services, interaction and annotations will rely on.
- In the current repository, the libraries under `src/` (core to io) are the kernel and its first supporting modules. Their order is enforced by `tools/check_layering.py`, which rejects any upward dependency, and is the mechanism the acyclic-graph rule extends to new modules. The Phase 0 debug viewer (TODO.md §8) is the beginning of the inspection application.

This proposal supplements these decisions and supersedes none.
