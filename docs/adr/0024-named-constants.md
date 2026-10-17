# ADR-0024: Named constants, no magic numbers

- **Status:** Accepted
- **Date:** 2026-10-12
- **Accepted:** 2026-10-17 by Onur Tuncer, on review; `tools/check_literals.py` passes over the tree and runs in CI
- **Owner:** Onur Tuncer
- **Phase:** 0
- **Related:** ADR-0001 (tolerances through the `Context`), ADR-0002 (predicate error bounds), ADR-0017 (naming), ADR-0021 (scalar-generic code), ADR-0022 (complexity limits and their lint)

## Context

Numeric literals in geometric code carry meaning that the digits do not show. `3.141592653589793` is π to the reader who recognises it and noise to the one who does not; `1e-5` is the model resolution in one place and a convergence threshold in another; `3.3306690738754706e-16` is Shewchuk's orientation bound :math:`(3 + 16u)u`, which nobody can verify from the decimal. A literal cannot be found by name, cannot be changed in one place, and in scalar-generic code (ADR-0021) cannot take the type of the scalar.

ADR-0001 already bans literal tolerances in `topo`, `boolean` and `mesh`, enforced by `tools/check_epsilon.py`. The first predicate tests and benchmarks of Phase 0 nevertheless arrived with π spelled out, because nothing said otherwise. C++20 provides `<numbers>` with `std::numbers::pi_v<T>` and the other common constants, so there is no longer a reason to type them.

## Options

1. **Review convention only.** Costs nothing; the literal that slipped through review is how this ADR came about.
2. **Named constants with a lint for the recognisable cases.** Mathematical constants come from `<numbers>`, model and algorithm constants are named once with their provenance, derived constants are written as the formula they come from; a script flags the digits of the well-known constants and long decimal literals in library code.
3. **Ban all numeric literals except 0 and 1.** Unworkable: dimensions, indices, loop counts, test inputs and expected values are literals by nature.

## Decision

Adopt option 2.

- **Mathematical constants** are taken from `<numbers>`: `std::numbers::pi`, `e`, `sqrt2`, `ln2`, `phi` and the rest, and `std::numbers::pi_v<T>` in code generic over the scalar. They are never spelled as decimals, in library code, tests, benchmarks or Python. Degree and radian conversions are written as `std::numbers::pi / 180` rather than `0.0174532925`.
- **Model and algorithm constants** are defined once, as a named `constexpr` (ADR-0017: PascalCase), with a comment stating where the value comes from: `Context::DefaultResolution` (ADR-0001), `Context::SizeBoxHalfWidth` (ADR-0001), `DefaultPageSize` (ADR-0005), `Splitter` (Dekker). Code refers to the name, never to the value.
- **Derived constants** are written as the formula from their source, not as the evaluated decimal: `(3.0 + 16.0 * Epsilon) * Epsilon` for Shewchuk's bound, `0x1p27 + 1.0` for the splitter. The compiler evaluates them; the reader can check them against the paper.
- **Powers of two** are hexadecimal floating literals, `0x1p-53`, never their decimal expansion.
- **Tolerances** follow ADR-0001 unchanged: through `Context` and `ToleranceOf`, never literals in `topo`, `boolean` or `mesh`.
- **Literals that stay literals:** 0, 1, 2 and the small integers of dimensions, indices, strides, loop counts and bit widths; test inputs and expected values; benchmark sizes, which are named at the top of the file. Exact decimal fractions such as `0.5` and `0.25` are fine where they mean "half" and "a quarter".
- **Exceptions** are marked on the line, `// STEREON-LINT-ALLOW(literal): <reason>`, like the epsilon lint's marker, and reviewed as code.
- **Enforcement.** `tools/check_literals.py` reads the tracked C++ sources (`src/`, `tests/`, `bench/`, `python/`) and fails on two patterns: a floating literal that begins with the digits of a `<numbers>` constant, its common multiples (2π, π/2, π/4) or the degree conversion factors, anywhere; and a decimal floating literal with ten or more significant digits in `src/` and `bench/`, where a value that precise is a constant that should carry a name and a formula. Tests are exempt from the second rule, since expected values are literals by nature. The lint job and CTest run the script with a self-test, as they do `check_epsilon.py`.

## Consequences

- Constants read as what they are and change in one place; generic code gets the right precision for its scalar type.
- The lint catches only what a regular expression can recognise: the famous constants and suspiciously precise decimals. A domain constant written as `0.02` is still a review matter.
- A few lines get longer (`std::numbers::pi` for `3.14159`). Accepted.
- Hexadecimal float literals are unfamiliar to some readers; the comment on the constant names the power.

## Verification

- `tools/check_literals.py --self-test` passes, with cases for π, e, √2 and the degree factor in several spellings, a ten-digit decimal in `src/` flagged and the same digits accepted in `tests/`, hexadecimal floats accepted, and the marker honoured.
- `tools/check_literals.py` reports no violation on the tree at adoption.
- The lint job and CTest run both.
