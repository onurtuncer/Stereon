# ADR-0002: Numeric robustness strategy

- **Status:** Proposed
- **Date:** 2026-09-28
- **Owner:** Onur Tuncer
- **Phase:** 0

## Context

Floating-point round-off makes geometric predicates (orientation, point-in-polygon, side-of-plane) give inconsistent answers for near-degenerate input. Inconsistent predicates are the most common root cause of invalid boolean results. We need robust topological decisions without paying for exact arithmetic everywhere.

## Options

1. **Plain floating point with epsilons.** Fast, and the source of most kernel bugs.
2. **Exact arithmetic everywhere** (GMP rationals, CGAL-style exact kernels). Robust but 10–100× slower, and it does not extend to NURBS evaluation.
3. **Exact predicates with floating-point filters, floating point for constructions.** Adaptive-precision predicates decide topology; constructions (intersection points, curves) stay in double precision within the tolerance model of ADR-0001.

## Decision

Adopt option 3.

- `stereon::robust` provides Shewchuk-style adaptive predicates: `Orient2D`, `Orient3D`, `InCircle`, `InSphere`, plus expansion arithmetic for custom predicates.
- An `Interval` type (outward-rounded, using directed rounding or error-free transformations) serves as the fast filter for curve/surface predicates; ambiguous intervals fall back to higher precision or to subdivision.
- Every topological decision in `topo/`, `boolean/` and `mesh/` goes through a predicate in `robust/`. Direct sign tests on computed doubles are banned in those libraries (clang-tidy rule).
- Constructions return a value plus an error bound where one is cheap to compute.

## Consequences

- Topology code becomes consistent by construction: the same query always returns the same answer.
- Predicates for curved geometry (point vs NURBS surface) remain approximate; they combine interval filters with the tolerance model rather than being exact.
- We must maintain a predicate test suite with adversarial near-degenerate cases.

## Verification

Pass Shewchuk's adversarial predicate tests and a 10⁶-case random near-degenerate suite with zero inconsistent results. Filtered predicates must run within 1.5× of naive floating point on non-degenerate input.
