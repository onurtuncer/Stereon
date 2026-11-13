# ADR-0052: Derivatives of meetings of space curves

- **Status:** Accepted
- **Date:** 2026-11-10
- **Accepted:** 2026-11-10 by Onur Tuncer: option 3 now, option 1 later as a separately named product.
- **Owner:** Onur Tuncer
- **Phase:** 2 (v0.2 "Intersections", intersection linearisation, ADR-0021 §9 stage 4)
- **Related:** ADR-0001, ADR-0010, ADR-0021, ADR-0029

## Context

ADR-0021 §9 asks every intersection handler to expose the linearisation of its results and a regularity classification. For a curve against a surface, a branch of two surfaces and a point where two plane curves meet, the result is a point that solves a square system and moves smoothly with the geometry wherever that system is regular.

Two space curves are different. Since #240 and #242 they meet where their closest approach is within the Context's resolution (ADR-0001): a meeting is the critical point of the distance between the curves, with the distance below 1e-5 mm by default. A generic change of either curve, of any size above the resolution, moves the curves apart and the meeting disappears; the primal `Intersect` then reports nothing. A meeting of space curves is therefore an event of codimension one, as a tangency is, not a point that moves. Differentiating it needs a convention, and the convention changes what the derivative means, so it is a design decision rather than an implementation detail.

## Options

1. **Track the closest approach.** Differentiate (s, t) at the critical point of the distance, (A(s) − B(t)) · A′(s) = 0 and (A(s) − B(t)) · B′(t) = 0, two equations in two unknowns, regular wherever the curves are not parallel. Always defined, consistent with how the kernel finds meetings, and gives the gap's derivative for free. But the derivative then describes the closest approach, not a meeting: past the resolution the primal result has no point to compare it with, and the reported position (on A, on B, or between) needs its own convention.
2. **Keep the meeting.** Treat the meeting as a constraint and differentiate only along design changes that keep the curves touching, projected onto that subspace; answer `Nondifferentiable` otherwise. Faithful to the primal result, but almost every change from a real design leaves the subspace, and the projection needs a further convention for which parameters give way.
3. **Answer that there is no derivative.** `Intersect`'s products for two space curves answer `Nondifferentiable` for a result with a meeting or an overlap, saying why, and `Valid` with no values for a proved separation, which a small enough change keeps. Honest and cheap; it gives nothing to a caller whose edges come from space curves meeting.

## Decision

Option 3 for `Intersect`, now; option 1 later, as its own product.

- `IntersectJvp(a, b, intersection, context)` and `IntersectVjp` for two `Curve` answer `Nondifferentiable` for a result with any point or overlap, with the reason, and carry no values. A proved separation is `Valid`, with no values forward and zeros in reverse. No binding type exists for space curves until a product needs one.
- The smooth tracking of option 1 becomes a separately named product, `ClosestApproach(a, b)`, with its own forward and reverse products and its own tracking convention in ADR-0021's terms. It is not hidden inside `Intersect`'s derivatives, so `Intersect`'s products always describe what `Intersect` returns.
- Points of two plane curves are not affected: in the plane a crossing persists under any small change and has the products of #257.

## Consequences

- `Intersect`'s derivatives stay tied to its primal result, as for every other handler: no product describes a point the primal call would not report.
- A design whose edges come from space curves meeting gets no sensitivity until `ClosestApproach` exists; the status says so rather than giving a number.
- `ClosestApproach` must choose how it reports the position (on A, on B, or the midpoint) and must be verified against a closest-approach solver, not against `Intersect`, whose result vanishes once the gap exceeds the resolution.
- The v0.2 ledger's linearisation row lists `ClosestApproach` as remaining work.

## Verification

- `tests/intersect/PlanarLinearizationTests.cpp`, *a meeting of two space curves is an event, a separation keeps*: two crossing lines answer `Nondifferentiable` with no values; two lines apart answer `Valid`, no values forward and a zero in reverse, and refuse a non-empty adjoint; a B-spline against itself, an overlap, answers `Nondifferentiable`.
- For `ClosestApproach`, when it is built: closed forms for two lines, second-order Taylor decay and the adjoint identity on random curves against an independent closest-approach solver, and `Nondifferentiable` for parallel curves.

## Addendum 2026-11-12: ClosestApproach design, owner decision

The owner accepted, on 2026-11-12, the design of option 1's product. It changes none of the decisions above: `Intersect`'s products for two space curves stay as decided.

1. **Result.** `ClosestApproach` reports both closest points `A(s)` and `B(t)`, the parameters `s` and `t`, the gap vector `d = A(s) − B(t)` and its length `|d|`. It reports no single "position": not on A, on B, or the midpoint.
2. **Domain ends.** A minimum at an end of a parameter's domain holds that parameter fixed, and its products differentiate the remaining equation. A minimum that switches between an end and the interior is `Nondifferentiable`.
3. **Ties.** A tie between candidate minima that the primal cannot decide is `Indeterminate`. A set of minimisers proved exactly is `Nondifferentiable`.
4. **Meetings.** When the curves actually meet (the distance is within the resolution), the whole result is `Nondifferentiable`. `SensitivityResult` carries a single status, so the parameters do not get a separate `Valid` status.
5. **Supported domain.** v0.2 certifies lines, circles and ellipses. B-spline pairs answer `NotSupported` until a subdivision proof of their minimum exists.
6. **Near events.** `NearEventSine` (2⁻²⁰) is reused, applied to the smaller eigenvalue of the Hessian normalised by the tangents' lengths.
7. **Tracking convention.** `ClosestApproach` declares `TrackingConvention::UniqueProjection` in the capability catalogue. No core API changes.

Implementation: `src/intersect/ClosestApproach.cppm` (a new `Stereon.Intersect:ClosestApproach` partition), `ClosestApproach.cpp` (the primal) and `ClosestApproachLinearization.cpp` (the products); `docs/sensitivity.rst`, *Closest approach of two space curves*.

Verification for this addendum: `tests/intersect/ClosestApproachTests.cpp`. It checks:

- closed forms for two skew lines, a line and a circle, and minima at a segment's end and at an arc's end;
- an independent oracle: dense sampling, then Newton's method in `long double`, with tolerances sized by `std::numeric_limits<long double>::epsilon()` (which equals `double`'s on some platforms);
- the metamorphic relations of exchange, rigid motion, reversal and scaling;
- second-order Taylor decay and the adjoint identity;
- every status, and the input refused.
