# ADR-0057: Inline storage for analytic circle seam splitting

- **Status:** Proposed
- **Date:** 2026-11-27
- **Phase:** 2
- **Related:** [ADR-0056, C++26 feature usage](0056-cpp26-feature-usage.md)

## Context

[SeamPieces](../../src/intersect/ElementaryPairs.cpp) partitions one circle parameter period wherever either of its two pcurves crosses a chart seam. Its cut list and returned interval list currently use dynamic storage. Both capacities follow from the helper's inputs.

## Proposed decision

Use `std::inplace_vector` for both lists inside this private intersection path. Name the capacities and document their derivation beside the declarations.

| List | Capacity | Derivation |
| --- | ---: | --- |
| Cut parameters | 4 | The two period endpoints plus at most one optional seam parameter from each of two calls to `SeamOf`. |
| Result intervals | 3 | Each successive pair of sorted cuts yields at most one interval. Define this capacity as the cut capacity minus one. |

Retain the existing ordering and strict gap test (`gap > 0x1p-40`). Apply the hot-tier capacity and inline-footprint checks proposed in ADR-0056. The bound ends at the return from `SeamPieces`. Subsequent chart clipping in `CircleSpans` can split intervals and requires its own capacity analysis.

## Alternatives and consequences

Reserving space in each vector would avoid growth reallocations but still require dynamically allocated storage. A fixed array plus a count would encode the same bound while requiring the caller or helper to track the active range. An inline vector carries its own size and supports the existing iteration and sorting operations.

The change removes heap allocation for these two lists and gives each a fixed inline footprint. If a pcurve can later contribute multiple seam parameters per period, revisit both capacities before extending the helper.

## Verification

Exercise zero, one and two distinct interior cuts, plus coincident cuts. Check the ordered interval endpoints. Include gaps at and immediately above the existing threshold to verify the strict comparison. Compare allocation counts on the same enclosing intersection fixture before and after the change. Report any timing change separately from allocation counts.
