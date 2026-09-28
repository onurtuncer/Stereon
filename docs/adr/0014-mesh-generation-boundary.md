# ADR-0014: Mesh generation boundary

- **Status:** Proposed (to be accepted before Phase 5)
- **Date:** 2026-09-28
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
- **Volume meshing:** optional adapters that hand the surface mesh and tags to Gmsh and Netgen, built as separate targets so their licences stay out of the core (ADR-0008).
- **Export:** STL, OBJ, glTF, and a tagged surface-mesh format suitable for CFD meshers.

## Consequences

- The kernel guarantees watertight, tagged surface meshes — the part general meshers handle worst when starting from CAD.
- Volume meshing quality depends on external tools; acceptable.
- Our own CDT must be robust (exact predicates, ADR-0002).

## Verification

Phase 5 exit gate: analysis-grade meshes of all corpus solids are watertight (every edge shared by exactly two triangles) and are accepted by Gmsh and snappyHexMesh without repair.
