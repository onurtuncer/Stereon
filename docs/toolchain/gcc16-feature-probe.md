# GCC 16 C++26 feature probe

Input for the ADR-0016 review. Each feature was compiled and run in a separate small probe program.

- **Date:** 2026-10-08
- **Toolchain:** GCC 16.2.0 (MSYS2 UCRT64, Rev4), CMake 4.4.4, Ninja 1.13.2, GDB 18.1
- **Platform:** Windows 11, `x86_64-w64-mingw32`
- **Base flags:** `-std=c++26 -O1`

## Summary

| Feature | Result | Flag needed | Feature-test macro |
| --- | --- | --- | --- |
| Contracts (P2900): `pre`, `post`, `contract_assert`, custom violation handler | ✅ Works | None. On by default with `-std=c++26`; `-fcontracts` is accepted | `__cpp_contracts 202502L`, `__cpp_lib_contracts 202502L` |
| Static reflection (P2996): `^^T`, splices, `nonstatic_data_members_of`, `identifier_of`, `define_static_array` | ✅ Works, including the member-walking serialiser pattern for ADR-0012 | `-freflection` (without it, `std::meta` is not declared) | `__cpp_impl_reflection 202603L`, `__cpp_lib_reflection 202603L`, `__cpp_lib_define_static 202506L` |
| Expansion statements (P1306): `template for` over tuple, init-list, constexpr array and reflected members | ✅ Works | None | `__cpp_expansion_statements 202506L` |
| Pack indexing (P2662): type packs `Ts...[I]` and value packs `as...[I]` | ✅ Works; the `VisitBuiltin` pattern for ADR-0003 compiles and dispatches correctly | None | `__cpp_pack_indexing 202311L` |
| Modules: named module, `import std;` | ✅ Works | `-fmodules`; `std` built from `bits/std.cc` with `-fsearch-include-path`; **`-lstdc++exp` on MinGW** to link `std::print`/`std::println` | `__cpp_modules 201810L`, `__cpp_lib_modules 202207L` |

### Library features ADR-0016 lists

| Feature | `__cpp_lib_*` | Consequence |
| --- | --- | --- |
| `std::execution` (senders) | not defined | **Keep the stdexec polyfill** |
| `std::simd` | not defined | **Keep the xsimd polyfill** |
| `std::linalg` | not defined | Not available; not needed in Phase 0 |
| `std::inplace_vector` | `202603L` | Use directly, no polyfill |
| `std::function_ref` | `202603L` | Use directly, no polyfill |
| `submdspan` | `202603L` | Use directly, no polyfill |
| `std::expected` | `202211L` | Use directly (ADR-0006) |

## Contract evaluation semantics

Selected with `-fcontract-evaluation-semantic=`. The default is `enforce`.

| Semantic | Behaviour observed on a precondition violation |
| --- | --- |
| `enforce` | Calls the custom `handle_contract_violation` (semantic reported as 3), then terminates |
| `observe` | Calls the handler for the precondition and the postcondition, then continues to the end of `main` |
| `quick_enforce` | Terminates immediately **without** calling the handler |
| `ignore` | Checks are not evaluated |

Other GCC options:

- `-fcontracts-client-check=[none|pre|all]`: check contracts at the call site for non-virtual functions.
- `-fcontracts-definition-check=[on|off]`: check contracts inside the callee (default `on`).
- `-fcontract-checks-outlined`, `-fcontract-disable-optimized-checks`, `-fcontracts-conservative-ipa`.

The default `handle_contract_violation` is a weak symbol in `libstdc++exp`, not in `libstdc++`, so any program built with `-fcontracts` and a checking semantic must link `-lstdc++exp` on every platform. Without it, Linux and macOS fail at link time with an undefined reference to `handle_contract_violation`.

On Windows, a terminating contract violation exits the process with code 127 and prints `terminate called without an active exception`. Test harnesses and death tests should expect that. CTest sees the termination as a crash (exception 0xC0000409) and fails the test before matching its output, so a test that expects a violation installs its own `handle_contract_violation` and exits normally, as `tests/probes/ContractProbe.cpp` does.

## Known compiler bugs

| Bug | Trigger | Workaround |
| --- | --- | --- |
| Internal compiler error (segfault) in GCC 16.2.0 | A defaulted **hidden-friend** comparison (`friend constexpr bool operator==(const T&, const T&) = default;`) in a type exported from a module, used in a translation unit that imports the module. Happens in plain modules and partitions, with or without `-freflection` and contracts. | Declare defaulted comparisons as members: `constexpr bool operator==(const T&) const = default;`. Defaulted member comparisons work. |
| Internal compiler error (segfault) in GCC 16.2.0 on stdexec | Any translation unit compiled with `-fmodules` that includes stdexec, even with no module code: `#include <stdexec/__detail/__query.hpp>` alone reproduces it (`__query.hpp:177`, `forwarding_query_t::operator()`). Seen with stdexec `nvhpc-25.09`, `nvhpc-26.05` and `main` (2026-10-09). CMake adds `-fmodules` to every source it scans for modules. Plain `-std=c++26` builds of stdexec work. A stdexec header unit builds but cannot be used (`no class template named '__f' in 'struct stdexec::just_t'`). | `Stereon::Exec` is a header, `<Stereon/Core/Exec.hpp>`, not a module partition. Files that include it are excluded from module scanning with `stereon_exec_sources()`, so they cannot import modules. |
| Internal compiler error C1001 in MSVC 19.50 on stdexec | stdexec senders used through a named module that includes stdexec in its global module fragment. stdexec included directly builds and runs. | Same as for GCC. |
| False `-Wshadow` warning in GCC 16.2.0 on expansion statements | `template for (constexpr std::meta::info e : ...)`: every unrolled iteration warns that `e` shadows the previous one. | Build reflected tables with an ordinary loop in a `consteval` function instead, as `ToString(ErrorCode)` does in `src/core/Error.cppm`. |
| `std::source_location` is TU-local in a module interface in GCC 16.2.0 | A module interface unit that includes `<source_location>` in its global module fragment and uses `std::source_location` in an exported declaration fails with "exposes TU-local entity 'struct std::source_location'": libstdc++ declares the type with internal linkage there. | `Stereon::SourceLocation` in `src/core/Log.cppm` is the kernel's stand-in, built from the same `__builtin_FILE()`, `__builtin_FUNCTION()` and `__builtin_LINE()` default arguments. |
| False `-Wnull-dereference` in GCC 16.2.0 on `std::vector<std::shared_ptr<T>>` copy assignment | At -O3 with `-Werror`, a defaulted copy assignment of a class holding a `std::vector<std::shared_ptr<T>>` reports "potential null pointer dereference" inside `std::__shared_count::operator=`, inlined through `std::copy` in `vector::operator=`. Seen in the release build of `tests/core/ArenaTests.cpp`; the same warning fires in stdexec's intrusive queue (see `stereon_exec_sources()`). | `PagedArena` (`src/core/Arena.cppm`) assigns by copy-and-swap, which copy-constructs the pointer vector instead of assigning over the old elements. |
| "contract condition is not constant" in constant evaluation in GCC 16.2.0 | A `static_assert` (or `STATIC_REQUIRE`) evaluates an exported constexpr function whose call chain reaches a hot-path precondition along certain paths (first seen next to an `if consteval`, then on the predicates' power-of-two rescaling path, `Expansion::PushBackNonZero`). The same calls evaluated at run time with contracts enforced pass, so no contract is violated. | Test those paths at run time (`tests/robust/DomainTests.cpp`); keep compile-time checks for the paths that GCC accepts. |
| "failed to read compiled module cluster: Bad file data" in GCC 16.2.0 | A partition of one module (`Stereon.Surface:Cylinder`) instantiates class and function templates imported from another module (`BasicCircle<T, 3>` and its `Eval`, `Deriv`, `Project` from `Stereon.Curve`) in its own templates and in a `static_assert`. The module builds, but a unit that imports it and then names a core template in a partial specialisation (`tests/support/Formatting.hpp`) fails to read the compiled module file. Reproducible from a clean build; bisected to that partition. | `BasicCylinder` computes its sections itself instead of building a `BasicCircle`. Until fixed, a geometry module does not build its functions out of another module's geometry templates; shared arithmetic belongs in `Stereon.Core`. |
| "failed to load binding '::sin@Stereon.Surface'" ("Bad file data") in GCC 16.2.0 on glibc | `Stereon.Surface:Cylinder` included `<cmath>` in its global module fragment and called `std::sin`, `std::cos` and `std::atan2`, while `Stereon.Surface` also imports `Stereon.Curve`, whose `Circle` partition includes `<cmath>` too. In the Linux optimised presets (`gcc16-release`, `gcc16-tsan`) a test that includes `<cmath>` after importing the module fails to read the module file. Debug, ASan, macOS and MinGW are unaffected. | `Stereon.Core:Trigonometry` provides `Sin`, `Cos` and `Atan2`, with the `double` overloads defined in `Trigonometry.cpp`, so no interface outside core refers to the C library's functions; geometry modules above core call those and do not include `<cmath>`. Confirmed with a Linux CI probe run. |
| Link error `undefined reference to __tu_has_violation` in GCC 16.2.0 | A unit imports a module whose **non-template inline** function carries a contract, and the unit parses no contract of its own. The check calls a per-unit helper that the importer never emits. | Compile one contract check in every importing unit, and build without optimisation: `Stereon::Detail::ContractAnchor()` in `<Stereon/Core/Contracts.hpp>`, which `stereon_anchor_contracts()` preincludes into Stereon's non-module targets, defines the helper at -O0. From -O1 up the importer inlines or clones the helper (`.constprop`) and never defines the plain name, so the anchor fails; `-fcontract-disable-optimized-checks -fcontract-checks-outlined` together restore it but warn "'optimize' attribute currently unsupported in modules" once per contract-bearing non-template function in module units, with no `-W` option to silence it, which `-Werror` CI cannot accept. Stereon therefore puts the contracts of inline module functions in the hot-path tier, which optimised presets compile out; `gcc16-observe` builds without optimisation; other importers include the header and build at -O0 or without hot contracts. |
| Contracts in **templates** imported from a module are silently ignored in GCC 16.2.0 | A unit instantiates a template with `pre`, `post` or `contract_assert` from an imported module: the check is compiled out whatever the unit's evaluation semantic, with or without the anchor, with run-time operands, at -O0 and at -O2. | None. Such checks fire only inside the module's own units. `stereon_contract_probe` (tests/probes) runs as a known-failing CTest that turns red when a compiler release fixes this; remove the anchor then. |

Report the bugs upstream and recheck them with each compiler release. A regression test should guard each workaround once the fix lands.

**Consequence for ADR-0007.** Until GCC's stdexec crash is fixed, code that builds sender graphs cannot live in module units. It has to sit in non-module source files behind plain function declarations, which module code can call through an ordinary header. This should be settled before Phase 1 algorithms are written.

## Findings for the ADR-0016 review

1. **Reflection is no longer partial.** GCC 16.2 reports the P2996 macros at `202603L`, and the serialiser probe walks struct members end to end. The ADR's Context section ("GCC 16 is the only compiler with reflection (partial)") can be updated.
2. **The release split in ADR-0016 needs a mechanism.** *Resolved: option 3, contract macros (see ADR-0016, "Mechanism for the release split").* The evaluation semantic is chosen per translation unit, not per contract, so "ignore for hot-path checks, enforce for cheap API-boundary checks" cannot be expressed with one global flag. Options to decide between:
   - Compile public API-boundary translation units with `enforce` or `quick_enforce` and internal ones with `ignore`. This is the mixed mode that P2900 allows.
   - Use `-fcontracts-client-check=pre`, so callers check preconditions while definitions are built with `-fcontracts-definition-check=off`.
   - Mark hot-path checks with a project macro that expands to nothing in release builds.
3. **Use `quick_enforce` for release boundary checks.** It is the cheapest terminating semantic, but it skips the violation handler, so release builds will not log diagnostics.
4. **Polyfills still needed:** stdexec for `std::execution`, xsimd for `std::simd`. The feature-test-macro wrappers (`Stereon::Exec`, `Stereon::Simd`) should switch on `__cpp_lib_senders` and `__cpp_lib_simd`. *Done: `Stereon::Simd` is the `Stereon.Core:Simd` partition; `Stereon::Exec` is a header because of the stdexec crashes above.*
5. **Build system notes for CMake:**
   - add `-freflection` globally
   - modules need `-fmodules` (CMake's C++26 module support adds this)
   - add `-lstdc++exp` on every platform: contracts need its default violation handler, and MinGW also needs it for `std::print`
   - `import std;` needs CMake's experimental `CMAKE_CXX_MODULE_STD` support, or building `bits/std.cc` ourselves

## Other compilers (non-blocking CI, 2026-10-08)

| Compiler | Builds and tests the Phase 0 skeleton | `-freflection` | `-fcontracts` | Notes |
| --- | --- | --- | --- | --- |
| Clang 22.1.8 (apt.llvm.org, libc++) | ✅ | Not accepted | Not accepted | Named modules work with CMake 3.31 and `clang-scan-deps-22` |
| MSVC 19.51 (windows-2025 runner) | ✅ | n/a | n/a | CMake 4.4 does not know C++26 for MSVC, so it builds with `CXX_STANDARD 23`, which selects MSVC's newest standard switch |

Neither can build code that uses reflection or contracts. Both stay non-blocking until they can (ADR-0016).

Until then, Stereon code keeps two fallbacks, both to be removed when these compilers catch up:

- **Contracts.** The boundary macros `STEREON_PRE`, `STEREON_POST` and `STEREON_ASSERT` (`<Stereon/Core/Contracts.hpp>`) exist only so that Clang and MSVC compile code with contracts. Once both accept `-fcontracts`, write `pre`, `post` and `contract_assert` directly. The `STEREON_HOT_*` macros stay, because C++26 has no per-contract labels to separate hot-path from boundary checks.
- **Contracts across module boundaries.** Every non-module target that imports a Stereon module is anchored (see the bug table); the contracts of inline module functions are hot-path contracts, which the optimised presets (`release`, `windows`, `bench`, `tsan`) compile out, and `gcc16-observe` builds without optimisation. A precondition on a template exported from `Stereon.Core` is enforced only for callers inside the module until GCC checks imported template contracts. A check that must fire in every importer at every optimisation level, such as the checked-handle report, is made from a non-inline function in a module implementation unit (`src/core/Handle.cpp`), whose contract is compiled in the module's own unit.
- **Reflection.** Code that uses reflection checks `__cpp_impl_reflection` and `__cpp_lib_reflection` and provides a reduced fallback; `ToString(ErrorCode)` returns `"ErrorCode"` without it.
