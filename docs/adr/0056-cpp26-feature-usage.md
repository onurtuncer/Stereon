# ADR-0056: Where C++26 features are used in kernel code

- **Status:** Proposed
- **Date:** 2026-11-14
- **Owner:** Onur Tuncer
- **Phase:** 0
- **Related:** [ADR-0007](0007-threading-and-determinism.md), [ADR-0012](0012-native-file-format.md), [ADR-0015](0015-ai-contribution-policy.md), [ADR-0016](0016-cpp26-adoption-and-contracts.md), [ADR-0019](0019-library-linking-packaging-and-abi.md), [ADR-0020](0020-cpp-c-abi-and-python-bindings.md), [ADR-0024](0024-named-constants.md), [ADR-0042](0042-cpp26-contract-validation-responsibilities.md), [ADR-0043](0043-optional-performance-instrumentation.md), [ADR-0046](0046-ownership-borrowing-and-pointer-lifetimes.md), [ADR-0054](0054-topological-event-detection.md)

This ADR complements ADR-0016, which decides *that* Stereon uses C++26 and how contracts behave. This ADR decides *where* each feature is used. It changes nothing in ADR-0007, ADR-0016 or ADR-0042 and supersedes none of them.

## Context

The kernel uses C++26 features unevenly, and the same pattern is written two ways in neighbouring libraries. ADR-0016 lists the features Stereon targets and the GCC 16 probe (`docs/toolchain/gcc16-feature-probe.md`) records which ones the toolchain ships, but neither says where a feature belongs in kernel code.

The case that forces the decision is `std::inplace_vector`. Curve projection keeps its Bernstein temporaries in an `inplace_vector` of 75 coefficients, and `docs/performance.rst` records the reasoning. The interval subdivision in `src/intersect` (`Bernstein.hpp`, `BezierArc.hpp`, `BezierPatch.hpp`) has the same shape and allocates on the heap at every halving and every evaluation. Pull request #296 removes part of that. Nothing told its author which of the two patterns was the rule.

Use on `main` at `ca74138`, counted by text search over `src`, `tests` and `bench`:

| Feature | Shipped by GCC 16.2 | Use today |
| --- | --- | --- |
| Contracts | Yes | About 100 macro uses in 20 files; policy in ADR-0016 and ADR-0042 |
| `std::inplace_vector` | Yes | 4 files: curve projection, ray classification, two analytic intersection handlers |
| Static reflection | Yes, with `-freflection` | `Serialize` and `Enum` in `core` |
| Expansion statements | Yes | 4 files, the per-degree kernels behind `DispatchFixed` |
| Pack indexing | Yes | 1 file |
| `std::simd` | No; xsimd behind `Stereon::Simd` | Curve batch evaluation and surface grid evaluation |
| `std::execution` | No; stdexec behind `Detail::ParallelFor` | `boolean`, `curve`, `surface` |
| `std::function_ref` | Yes | None |
| `submdspan` and `mdspan` | Yes | None |
| `std::linalg` | No | None |

Three constraints shape any rule:

- **Determinism.** ADR-0007 promises bit-identical results for any scheduler and thread count (G1). A rewrite that reorders floating-point operations breaks it.
- **Review cost.** `intersect`, `boolean` and `blend` are Tier C under ADR-0015, so every change there is reviewed line by line.
- **Toolchain.** GCC 16.2 cannot compile stdexec in module units and drops contracts in imported templates (ADR-0016).

## Options

1. **No usage policy; leave it to review.** No rule to maintain. The cost is what exists today: one pattern in `curve`, another in `intersect`, and reviewers with no criterion to cite.
2. **Prefer the C++26 facility wherever it fits.** Simple to state. It invites stack growth (a 26 x 26 pole net held inline is about 43 KB), reordered arithmetic, and churn in Tier C code with no measured gain.
3. **Per-feature rules.** Each feature gets named places where it is used, places where it is not, the route (standard name or Stereon wrapper), and the evidence a change must bring. More text, but each rule is checkable.

## Decision

Adopt option 3. Three rules apply to every feature, then each feature has its own.

### Rules for every feature

- **Route.** A feature libstdc++ ships is used under its standard name. A feature it does not ship is used only through the Stereon wrapper (`Stereon::Simd`, `Detail::ParallelFor`). A feature with neither is not used.
- **No result changes.** A change that adopts a feature keeps the floating-point operations and their order. Existing tests that compare with `==` must still pass unchanged (ADR-0007, G1).
- **Evidence for speed.** A change made for performance adds or names a `HotPath` benchmark and a line in `docs/performance.rst`. It records no figure until the reference runner has produced one (ADR-0043).

### `std::inplace_vector`

Use it for a temporary whose largest size follows from a kernel bound.

- **Capacity is a named constant with its derivation.** For example `BSplineMaxDegree + 1` poles for one Bezier segment, or 4 tracks for two skew cylinders. No literal capacities (ADR-0024).
- **Capacity is the proved maximum, not a typical size.** There is no heap fallback. A small-buffer type that spills to the heap needs its own ADR.
- **The bound is a contract.** Code that fills the container from a run-time size states `STEREON_HOT_ASSERT(size <= Capacity)` first. If the size can come from outside the kernel, checked code returns a `KernelError` before that point (ADR-0016). No code relies on `std::bad_alloc` from an overflow.
- **Footprint limit.** One object holds at most 4 KiB inline. Larger data stays on the heap and is split in place with a stride.
- **Not for stored state.** Members of objects that outlive the call stay `std::vector`, so those objects stay compact. Curve projection's stored polynomials are the model.
- **Not in public interfaces.** No `inplace_vector` in exported signatures, the C ABI or the bindings (ADR-0019, ADR-0020).
- **Recursion.** In recursive subdivision, the inline footprint of one frame times the depth limit is stated in a comment beside the capacity.

### Contracts

ADR-0016 and ADR-0042 govern contracts. This ADR adds one obligation: every bound that an `inplace_vector` capacity, a `Simd` lane count or a `submdspan` extent relies on is stated as a hot-tier contract at the point of use.

### `std::simd`

- **Only through `Stereon::Simd`.** No `std::simd` or xsimd names outside `Stereon.Core:Simd`.
- **Only where values are computed, not decided.** Evaluators and batch or grid kernels (ADR-0007). Not in predicates, sign tests or interval enclosures, until a written argument shows the lanes round as the scalar code does.
- **Same arithmetic per lane.** A vector kernel performs the scalar kernel's operations in the scalar order, and the tail uses the same kernel. No reduction whose order depends on the lane width.

### `std::execution`

- **Only through `Detail::ParallelFor` on the `Context`'s pool.** No sender code in module units, and no operation creates threads (ADR-0007, ADR-0016).
- **Parallelise across independent items.** Face pairs, patch pairs, rays. Not inside one subdivision recursion, unless a benchmark shows a gain and the merge order is fixed by input index.

### Reflection, expansion statements, pack indexing

- **For mechanical code only.** Serialisers (ADR-0012), enumerator names, diagnostic structure (ADR-0042), bindings (ADR-0020), and per-size kernels through `DispatchFixed`.
- **Compile time only.** No reflection result is consulted on a hot path at run time.
- **Meaning stays explicit.** Semantic rules and stable diagnostic identifiers are written by hand; reflection generates the structure around them (ADR-0042).

### `std::function_ref`

- **For a callable that is passed down and not stored.** Loop bodies, visitors, kernels. It is a borrow under ADR-0046.
- **`std::function` only for a callable that is stored.** New parameters of type `const std::function&` are not added.

### `mdspan` and `submdspan`

- **For views of pole nets and sample grids.** They replace vectors of vectors and hand-written stride arithmetic in new code. Layout is row-major and stated where the view is made.
- **Non-owning and call-scoped.** A view is not stored beyond the call that made it (ADR-0046).

### `std::linalg`

Not used. libstdc++ does not ship it and Stereon adds no polyfill. Small fixed-size algebra stays in Stereon's own types. Revisit when `__cpp_lib_linalg` is defined, and then only for dense solves whose results are not compared bitwise.

### Features not listed

A C++26 feature not named above is used in kernel code only after this ADR is amended with its rule.

## Consequences

**Easier**

- Review has criteria to cite: route, bound, footprint, evidence.
- Subdivision loops in `intersect` can be made allocation-free the way projection already is. ADR-0054 names adaptive subdivision and interval arithmetic among the techniques for event detection.
- Moving from xsimd and stdexec to the standard library stays a change in one wrapper each.

**Harder**

- Capacity constants depend on `BSplineMaxDegree`. Raising it grows every dependent stack frame, so the footprint comments must be rechecked.
- Converting existing Tier C code needs line-by-line review for a change that alters no result (ADR-0015).
- `std::function` parameters and vectors of vectors that exist today become known debt. This ADR does not require converting them; it stops new ones.

**To build**

- A lint check, `tools/check_cpp26_usage.py`, for the rules a text search can decide: direct `std::simd`, xsimd or stdexec names outside the wrappers, literal `inplace_vector` capacities, new `const std::function&` parameters.
- An allocation-count harness for tests. The repository has none, so "this evaluation does not allocate" cannot be asserted today. Pull request #296 had to state that gap.
- `HotPath` benchmarks for the subdivision paths; #296 adds `HotPath/Intersect/SplinePair`.

## Verification

| Criterion | Pass condition | Status |
| --- | --- | --- |
| Lint check | `check_cpp26_usage.py` reports no violation in the lint job | Not built |
| No result changes | A conversion passes the existing suite, the OCCT differential tests and the corpus with no test edited | Met for #296 on Linux and macOS, including asan and tsan |
| No allocation in a converted evaluation | An allocation-count test reports zero for `EvalH` and for curve projection's Bernstein evaluation | Blocked: no harness |
| Speed | The named `HotPath` benchmark is not slower on the reference runner, by the method of ADR-0043 | No figures yet |
| Footprint | No `inplace_vector` type in `src` exceeds 4 KiB, checked by a `static_assert` beside each capacity | Not built |

Pull request #296 is the first change made under these rules and is still a draft. Its one failing job, Windows `gcc16-debug`, fails a randomly seeded plane and sphere property test in code the change does not touch.

### Open questions

- **The 4 KiB limit.** It is a proposal sized so that 26 interval poles (about 1.7 KB) and 75 Bernstein coefficients (600 bytes) pass and a pole net does not. It needs a number from the benchmark, or the owner's judgement.
- **Converting `Coefficients`, `Arc::Poles` and `Segment::Poles`.** These are types shared by eight handler units. Decide after `HotPath/Intersect/SplinePair` has figures.
- **`mdspan` in `BezierPatch.hpp`.** Whether the pole net moves from a vector of vectors to one buffer with a view now, or with the ADR-0054 work.
- **SIMD in interval code.** Whether batched box and normal-cone tests are worth a rounding argument, or stay scalar.
- **Reflection and the file format.** ADR-0012 generates serialisers by reflection. Whether a reflected member order may define a persistent layout without an explicit check is not settled here.
