# ADR-0014: Mesh generation boundary

- **Status:** Proposed (to be accepted before Phase 5)
- **Date:** 2026-09-28
- **Revised:** 2026-10-10 (tag provenance; end-to-end example)
- **Owner:** Onur Tuncer
- **Phase:** 5

## Context

Two meshing needs exist: display tessellation for viewers, and analysis-grade watertight surface meshes as input to CFD/FEA meshers. Volume meshing is a large field with mature open-source tools (Gmsh, Netgen, TetGen, snappyHexMesh). We must decide what the kernel owns.

## Options

1. **Surface tessellation only; no analysis meshing.** Smallest scope; analysis users must repair meshes themselves.
2. **Own surface meshing (display and analysis-grade); volume meshing via external tools.**
3. **Full in-house surface and volume meshing.** Maximum control; years of extra work outside the kernel's core competence.

## Decision

Adopt option 2.

- **Display tessellation:** chordal and angular tolerance controls, parallel per face, constrained Delaunay triangulation (CDT) in parameter space.
- **Analysis-grade surface mesh:** shared edge discretisation between adjacent faces guarantees watertightness; size fields from curvature and user controls; quality targets (minimum angle, aspect ratio); face and edge tags preserved for boundary conditions.
- **Tags come from topology, not from the mesher.** Boundary-condition tags are sparse attributes on faces and edges (ADR-0004) that survive booleans through the propagation policy of ADR-0011. The mesher copies each triangle's tag from its face and never invents or merges tags; a face without a tag yields untagged triangles, which the export reports. The mesh also carries the face ID of every triangle, so a result can be traced back to the B-rep.
- **The reference example.** The first end-to-end use case the kernel must pass, from Phase 4 onwards, is: an engineering duct body with `inlet`, `outlet` and `wall` tags, cut by a cylinder (boolean, ADR-0013), checked at the `Solid` level (ADR-0018), meshed watertight with tags, and accepted by an external volume mesher. It exercises tolerances, intersection completeness, classification, history, validity and tagging in one run.
- **Volume meshing:** optional adapters that hand the surface mesh and tags to Gmsh and Netgen, built as separate targets so their licences stay out of the core (ADR-0008).
- **Export:** STL, OBJ, glTF, and a tagged surface-mesh format suitable for CFD meshers.

## Consequences

- The kernel guarantees watertight, tagged surface meshes — the part general meshers handle worst when starting from CAD.
- Volume meshing quality depends on external tools; acceptable.
- Our own CDT must be robust (exact predicates, ADR-0002).

## Verification

Phase 5 exit gate: analysis-grade meshes of all corpus solids are watertight (every edge shared by exactly two triangles) and are accepted by Gmsh and snappyHexMesh without repair. The duct example produces a mesh in which every triangle carries exactly one of `inlet`, `outlet`, `wall`, the hole wall is `wall`, the inlet and outlet patch areas match the B-rep face areas within the chordal tolerance, and Gmsh produces a volume mesh with those three physical groups.
