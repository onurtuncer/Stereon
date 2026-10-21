# NOT YET REPRODUCED. Do not file as is.
# Summary (if it reproduces): [modules] std::source_location in exported decl: "exposes TU-local entity"
Component: c++ or libstdc++ (undecided)   Keywords: rejects-valid
Known to fail: 16.2.0 (seen in the Stereon build, compiler and target not recorded; probably MSYS2 UCRT64)
Not reproduced on: Compiler Explorer GCC 16.2.0 and trunk (x86_64 Linux), MinGW gcc 16.2.0 (CE)
Blocks: c++-modules (PR103524)

## Description
Reported in Stereon (docs/toolchain/gcc16-feature-probe.md): a module
interface that includes <source_location> in its global module fragment and
uses std::source_location in an exported declaration was rejected with
"exposes TU-local entity 'struct std::source_location'". The original
failing code was never committed; the repository only has the workaround
(Stereon::SourceLocation in src/core/Log.cppm, commit 9afc897, PR #25).

## Testcase (candidate; builds and runs correctly on CE)
m.cppm:
```c++
module;
#include <source_location>
export module m;
export struct Record {
  std::source_location loc;
};
export inline Record make(std::source_location loc = std::source_location::current()) {
  return {loc};
}
```

main.cpp:
```c++
import m;
int main() { return make().loc.line() == 2 ? 0 : 1; }
```

## Command lines
```
g++ -std=c++26 -fmodules -Wall -Wextra -c m.cppm -o m.o
g++ -std=c++26 -fmodules -Wall -Wextra -c main.cpp -o main.o
g++ m.o main.o -o a.out
./a.out
```
<gcc -v output here>

## What was tried (all compile cleanly and run correctly)
On g162 unless noted, -std=c++26 -fmodules -Wall -Wextra:
1. Exported inline function with a `std::source_location = current()` default argument (also on gsnapshot).
2. Exported struct with a std::source_location member plus 1.
3. 2 with -freflection -fcontracts -fcontract-evaluation-semantic=enforce.
4. Exported class with a non-inline member, an inline member and a member
   template, all taking `std::source_location = current()`; a struct
   with an NSDMI `= current()`; an exported constexpr variable, a constexpr
   function returning current(), an inline variable.
5. 4 with <format>, <meta> and <functional> also in the GMF, -freflection
   -fcontracts -O2.
6. `import std;` (--compile-std-module) instead of the GMF include, with
   a static constexpr Current(std::source_location = current()) member.
7. A Stereon-like partition (:log) re-exported by the primary interface:
   enum, LogRecord {level, string_view, source_location = current()},
   std::function sink, a Ctx class with a variadic std::format_string Log
   template and Emit(..., std::source_location = current()), with -freflection.
8. Exported `inline constexpr` and `constexpr auto` source_location
   variables, a lambda with a source_location default argument, a variable
   template, `constexpr const char* = current().file_name()`, a function-local
   `static constexpr auto = current()`.
9. 2 and 7 on MinGW gcc 16.2.0 (mingw64_ucrt_gcc_1620; compiled and linked, not run).
10. GMF include of <source_location> together with `import std;` in the
    interface and the importer.

## Notes
- libstdc++ 16.2 declares `struct source_location` as an ordinary
  namespace-scope class in std (it has a private nested `struct __impl`, a
  consteval `current(__builtin_ret_type = __builtin_source_location())`, and
  befriends `std::contracts::contract_violation`). Nothing in it has internal
  linkage, so the feature-probe explanation ("libstdc++ declares the type with
  internal linkage there") is not supported by the header.
- If it does reproduce locally, -fdiagnostics-color=never output and the
  exact declaration it names are needed to decide c++ vs libstdc++.
- Unrelated finding: case 8, `file != nullptr` in the importer, gives a
  garbled warning: "the address of ''string_cst' not supported by
  dump_decl<declaration error>' will never be NULL [-Waddress]". It is a
  separate, minor diagnostic bug.

## Duplicate check   (remove before filing)
WebSearch (Bugzilla is behind an Anubis check): "gcc bugzilla modules
source_location TU-local" found PR104602 and PR105128 (source_location,
not modules) and a gcc-help thread from Feb 2025 on TU-local alias exposure in
{fmt}. No duplicate.
Final Bugzilla search to be done by hand: `source_location TU-local`,
`exposes TU-local entity source_location`, `modules source_location`.
