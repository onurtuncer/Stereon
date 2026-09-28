# ADR-0016: C++26 adoption and contracts policy

- **Status:** Proposed
- **Date:** 2026-09-28
- **Owner:** Onur Tuncer
- **Phase:** 0

## Context

Stereon targets C++26 for static reflection, contracts, `std::execution`, `std::simd`, `std::inplace_vector`, `std::function_ref`, `submdspan`, pack indexing and expansion statements. As of September 2026, GCC 16 is the only compiler with reflection (partial) and contracts; Clang and MSVC lag, and no standard library ships `std::execution` or `std::linalg`. We must decide which features are mandatory and how contracts behave in each build.

## Options

1. **Stay on C++23** until all three major compilers support C++26 fully. Portable; loses reflection and contracts for years.
2. **Require full C++26 on all compilers.** Not buildable today.
3. **C++26 with a reference compiler, polyfills for missing library features, and non-blocking CI for other compilers.**

## Decision

Adopt option 3.

- **Reference compiler:** GCC 16+. CI on GCC is blocking; Clang and MSVC run as non-blocking until they support reflection and contracts, then become blocking.
- **Mandatory language features:** contracts, static reflection, expansion statements, pack indexing, C++20 modules.
- **Library features with polyfills** behind thin `stereon::` wrappers: `std::execution` → NVIDIA stdexec; `std::simd` → xsimd where libstdc++ support is incomplete. Wrappers switch to the standard version via feature-test macros.
- **Contracts.** Every public geometric operation states preconditions (parameter in domain, non-degenerate input, knot vectors non-decreasing) and key postconditions (finite results, valid output). Evaluation semantic per build:
  - `debug` and fuzzing: **enforce** (terminate with diagnostics);
  - corpus and CI runs: **observe** (log and continue, then fail the run);
  - `release`: **ignore** for hot-path checks, **enforce** for cheap API-boundary checks.
- Contract violations are programming errors and never become `KernelError`s (ADR-0006).

## Consequences

- Reflection removes serialisation and binding boilerplate (ADR-0012); contracts make interfaces self-documenting and catch misuse early.
- The project is tied to GCC for now, which limits who can build it; accepted as a deliberate trade-off.
- Polyfill wrappers add a thin maintenance layer until the standard libraries catch up.

## Verification

CI confirms the GCC build with contracts enforced passes all tests; feature-test macros correctly select standard vs polyfill implementations; the contract-observe corpus run reports zero violations before each release.
