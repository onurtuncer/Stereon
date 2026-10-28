# ADR-0045: CppAD as the automatic-differentiation backend

- **Status:** Proposed
- **Date:** 2026-10-28
- **Owner:** Onur Tuncer
- **Phase:** 0
- **Related:** ADR-0008 (licence), ADR-0016 (C++26), ADR-0019 (packaging), ADR-0020 (bindings), ADR-0021 (shape derivatives)

## Context

ADR-0021 makes Stereon's numerical geometry scalar-generic. Poles, weights,
radii and transform parameters may be any `Scalar`, while knots, parameters,
topology and robust predicates stay `double`. ADR-0021 selects no AD
backend. Its tests use a forward-mode dual number written for them
(`tests/support/Dual.hpp`), and the first derivative products (`EvalJvp`,
`EvalVjp`) are hand-derived. A kernel audit asked for evidence that the
generic evaluators work with a production AD library: taped reverse mode,
Jacobians, and composed transformations. The owner chose CppAD.

CppAD (COIN-OR, <https://github.com/coin-or/CppAD>) is a mature
operator-overloading AD library. It offers forward and reverse mode of any
order, sparse Jacobians and Hessians, and taping with retaping checks. It is
distributed under EPL-2.0, with GPL-2.0-or-later as a secondary licence.

EPL-2.0 is not on ADR-0008's allow-list for the libraries. ADR-0008 admits
copyleft dependencies "only as optional, separately built adapters", and
ADR-0021 asks for AD-library-specific behaviour behind internal adapters, as
a private dependency (ADR-0019).

## Options

1. **CppAD as an optional adapter.** It is mature, fast and complete, and it
   works with the existing scalar-generic code unchanged, because
   `CppAD::AD<double>` satisfies `Stereon::Scalar`. Its licence keeps it out
   of the libraries.
2. **A permissively licensed backend (autodiff, Enzyme, CoDiPack).** These
   could be linked into the libraries. autodiff (MIT) records no tape that
   can be replayed at new inputs, Enzyme is an LLVM plugin and Stereon is
   built with GCC, and CoDiPack is GPL-3.0, so not permissive either. The
   owner chose CppAD.
3. **No backend; hand-derived products only.** Exact and dependency-free,
   but every operation needs a derivation, and no taped reverse mode is
   available to applications.

## Decision

Adopt option 1. CppAD is Stereon's production AD backend, used as an
optional adapter under ADR-0008.

- **Pinned and optional.** CppAD release 20260000.0 is fetched at a pinned
  commit, only when `STEREON_WITH_CPPAD` is on (off by default). It is
  configured and installed into the build tree by its own CMake, which
  generates `cppad/configure.hpp` and builds `cppad_lib`, the small static
  library its headers call into.
- **Never in a library.** No target under `src/` links CppAD, and no source
  under `src/` includes it. `tools/check_dependencies.py` enforces both, and
  checks the fetched licence text against the declaration. The libraries
  stay scalar-generic and backend-neutral. An application that wants CppAD
  instantiates the evaluators with `CppAD::AD<double>` in its own code,
  under CppAD's licence.
- **What the kernel promises.** The scalar-generic evaluators compile with
  `CppAD::AD<double>`, and the tapes they record are valid for every value
  of the active inputs: their branches depend on knots and parameters only,
  never on poles, weights or transform parameters (ADR-0021, section 2).
  Taping an operation whose branches do depend on active values, such as a
  projection or an intersection, is outside this promise. Those are
  differentiated through their equations (ADR-0021, section 6).
- **Products stay native.** The kernel's own products (`EvalJvp`, `EvalVjp`
  and those to come) keep their hand-derived, `double` implementations and
  their `SensitivityResult` contract. CppAD verifies them, and serves
  applications that tape larger computations around them.
- **Bindings.** Python and the C ABI expose derivative products on numeric
  buffers (ADR-0020), never CppAD types.

## Consequences

- Applications get taped forward and reverse mode, Jacobians and Hessians
  of Stereon's generic evaluators without a second implementation.
- Each new scalar-generic path gains a CppAD integration test, so a change
  that branches on an active value, which would invalidate a tape, is
  caught.
- CppAD's licence obligations stay with the applications that link it. A
  binary distribution of Stereon's libraries carries no CppAD code.
- One more fetched dependency, built at configure time in the jobs that
  enable it.

## Verification

`tests/ad/CppADTests.cpp`, built with `STEREON_WITH_CPPAD=ON` and run in CI:

1. `CppAD::AD<double>` satisfies `Stereon::Scalar` (a `static_assert`).
2. The Jacobian CppAD tapes of B-spline curve points, in poles and weights,
   equals the closed form `R_i` and `N_i (P_i − C) / W` on random polynomial
   and rational curves of degree 1–5, to 1e-12 relative.
3. CppAD's forward mode equals forward-mode duals pushed through the same
   evaluator, and its reverse mode is the transpose of its forward mode.
4. One tape reproduces `Eval` for new poles and weights, and records no
   comparison of active values (`compare_change_number() == 0`).
5. A point moved by a composed rotation and translation has the closed-form
   derivative `k × R(x − c)` in the angle and the identity in the offset.
6. The Jacobian of a rational B-spline surface point in its poles is the
   tensor-product basis `N_i M_j w_ij / W`.
7. The kernel's own `EvalJvp` and `EvalVjp` equal CppAD's forward and
   reverse modes on the same curves, seeds and adjoints, to 1e-11 relative.
8. `tools/check_dependencies.py` fails when a library links or includes
   CppAD, and when the fetched licence is not the declared one.
