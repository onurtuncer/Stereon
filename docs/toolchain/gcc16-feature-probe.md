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

On Windows, a terminating contract violation exits the process with code 127 and prints `terminate called without an active exception`. Test harnesses and death tests should expect that.

## Known compiler bugs

| Bug | Trigger | Workaround |
| --- | --- | --- |
| Internal compiler error (segfault) in GCC 16.2.0 | A defaulted **hidden-friend** comparison (`friend constexpr bool operator==(const T&, const T&) = default;`) in a type exported from a module, used in a translation unit that imports the module. Happens in plain modules and partitions, with or without `-freflection` and contracts. | Declare defaulted comparisons as members: `constexpr bool operator==(const T&) const = default;`. Defaulted member comparisons work. |
| Internal compiler error (segfault) in GCC 16.2.0 on stdexec | Any translation unit compiled with `-fmodules` that includes stdexec, even with no module code: `#include <stdexec/__detail/__query.hpp>` alone reproduces it (`__query.hpp:177`, `forwarding_query_t::operator()`). Seen with stdexec `nvhpc-25.09`, `nvhpc-26.05` and `main` (2026-10-09). CMake adds `-fmodules` to every source it scans for modules. Plain `-std=c++26` builds of stdexec work. A stdexec header unit builds but cannot be used (`no class template named '__f' in 'struct stdexec::just_t'`). | `Stereon::Exec` is a header, `<Stereon/Core/Exec.hpp>`, not a module partition. Files that include it are excluded from module scanning with `stereon_exec_sources()`, so they cannot import modules. |
| Internal compiler error C1001 in MSVC 19.50 on stdexec | stdexec senders used through a named module that includes stdexec in its global module fragment. stdexec included directly builds and runs. | Same as for GCC. |

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
   - add `-lstdc++exp` on MinGW
   - `import std;` needs CMake's experimental `CMAKE_CXX_MODULE_STD` support, or building `bits/std.cc` ourselves

## Other compilers (non-blocking CI, 2026-10-08)

| Compiler | Builds and tests the Phase 0 skeleton | `-freflection` | `-fcontracts` | Notes |
| --- | --- | --- | --- | --- |
| Clang 22.1.8 (apt.llvm.org, libc++) | ✅ | Not accepted | Not accepted | Named modules work with CMake 3.31 and `clang-scan-deps-22` |
| MSVC 19.51 (windows-2025 runner) | ✅ | n/a | n/a | CMake 4.4 does not know C++26 for MSVC, so it builds with `CXX_STANDARD 23`, which selects MSVC's newest standard switch |

Neither can build code that uses reflection or contracts. Both stay non-blocking until they can (ADR-0016).
