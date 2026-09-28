# ADR-0009: Parameter-space curves (pcurves)

- **Status:** Proposed (to be accepted before Phase 1 topology work)
- **Date:** 2026-09-28
- **Owner:** Onur Tuncer
- **Phase:** 1

## Context

A face is a region of a surface bounded by edges. Trimming, point classification, tessellation and face splitting in booleans all work most reliably in the surface's (u, v) parameter space. Each edge therefore needs a 2D curve (pcurve) on every face it bounds, and the pcurves must agree with the 3D curve within the edge tolerance.

## Options

1. **Pcurves computed on demand** by projecting 3D curves. Less storage; projection is slow and can fail near singularities (cone apex, sphere poles).
2. **Pcurves mandatory on every coedge,** kept consistent with the 3D curve by construction.
3. **Pcurves only, 3D curves derived.** Compact; 3D curves from two faces rarely agree exactly.

## Decision

Adopt option 2.

- Every coedge references a pcurve on its face's surface. Every edge keeps its 3D curve.
- Operations that create edges (intersection, sweep, boolean splitting) produce the 3D curve and all pcurves together, with an error bound.
- The validity checker verifies that each pcurve mapped through its surface lies within the edge tolerance of the 3D curve.
- Periodic surfaces use explicit seam edges; pcurves never wrap silently across a period.
- Singular points (poles, apexes) are represented by degenerate edges with a pcurve and no 3D extent.

## Consequences

- Tessellation, trimming and classification work in 2D without projection.
- Storage roughly doubles for edge geometry; acceptable.
- Every edge-creating algorithm must output pcurves, which raises the bar for surface–surface intersection (ADR-0010).

## Verification

The validity checker passes on all Phase 3 constructed solids and on imported STEP files after healing; pcurve/3D deviation stays within edge tolerance at 1 000 samples per edge.
