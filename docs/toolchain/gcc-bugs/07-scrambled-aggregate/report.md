# Summary: [modules] wrong code: importer passes empty-class argument in a stack slot (x86-64)
Component: c++   Keywords: wrong-code, ABI
Known to fail: 12.2.0, 13.3.0, 14.3.0, 15.2.0, 16.1.0, 16.2.0, 17.0 (trunk 20261021)   Known to work: none
Blocks: c++-modules (PR103524)
Target: x86_64-pc-linux-gnu

## Description
On x86-64 (SysV psABI, -fabi-version >= 12), an argument of empty class type
occupies no stack slot. When a function is declared in a module and called
from a unit that imports the module, the importer still reserves a stack slot
for an empty-class argument that comes before an argument passed in memory.
The callee then reads every later memory argument 8 bytes too low. The
testcase below prints garbage instead of 123. The same code without modules
prints 123.

## Testcase
m.cppm:
```c++
export module m;
export struct Empty {};
export struct Big { long a, b, c; };
export struct S {
  static long f(Empty, Big b) { return b.a * 100 + b.b * 10 + b.c; }
};
```

main.cpp:
```c++
#include <cstdio>
import m;
int main() {
  long r = S::f(Empty{}, Big{1, 2, 3});
  std::printf("%ld\n", r);
  return r == 123 ? 0 : 1;
}
```

## Command lines
```
g++ -std=c++26 -fmodules -Wall -Wextra -c m.cppm -o m.o
g++ -std=c++26 -fmodules -Wall -Wextra -c main.cpp -o main.o
g++ m.o main.o -o a.out
./a.out
```
<gcc -v output here>

## Compiler output
No diagnostics. Output at -O0 (the value varies from run to run):
```
14073217128088812
```
Output at -O2:
```
12
```
Expected: `123`.

## Notes
- The callee sees `b` as {garbage, 1, 2}: shifted by one 8-byte slot.
- The importer is the side that gets it wrong. If the declarations are in an
  `export extern "C++" { }` block, an importer calling `T::h(Empty, Big)`
  defined in a non-module TU gets 112. A non-module TU calling `T::k`
  defined in the module interface gets 123.
- Fails: static member functions, non-static member functions, constructors
  (inline or defined out of line in the interface), and non-inline
  namespace-scope functions defined in the interface and called from the
  importer.
- Works: an inline namespace-scope function (caller and callee are both
  emitted in the importer and agree), an empty argument after the memory
  argument (`f(Big, Empty)`), and an empty argument before register-passed
  scalars (`f(Empty, long, long, long)`).
- It also fails within a single importer for inline member functions, so
  the importer's call site and its copy of the callee body disagree too.
- -O0 and -O2 both fail, and `-std=c++20 -fmodules-ts` fails the same way on
  12, 13, 14 and 15. It is not a regression.
- Seen originally in a larger project: a {uint32_t, std::string_view,
  uint32_t} aggregate passed by value to a constructor whose first parameter
  was an empty tag class (`std::make_shared<T>(Key{}, tables, cert, id)` in a
  module implementation unit) was stored with its fields shifted by 8 bytes.
  Only the Linux build failed. Windows x64 passes large aggregates by
  reference, so MinGW is unaffected; the macOS build was unaffected too.
- Workaround: do not pass empty tag classes by value ahead of memory-passed
  arguments to functions from a module (put them last, or pass them by
  reference).

## Duplicate check   (remove before filing)
Bugzilla itself is behind an Anubis bot check, so it was searched through
WebSearch only:
- "gcc bugzilla modules empty class parameter passing wrong-code TYPE_EMPTY_P"
  (gcc.gnu.org): only PR90164 (-Wabi wording) and the GCC 8/10 release-note
  ABI changes. Not a duplicate.
- "gcc modules wrong code empty struct argument imported member function ABI":
  PR84502 (GCC 8 empty class template, -fabi-version=12; no modules) and
  PR89627 (RISC-V flattening). Not duplicates.
- "modules empty class argument passed by value wrong-code import -fmodules
  stack argument shifted" (gcc.gnu.org): nothing relevant.
No likely duplicate found. The cause is probably that the empty-record flag
(TYPE_EMPTY_P / the target's empty_record_p decision) is not streamed for
imported classes, or that the importer's function types do not carry it.
Final Bugzilla search to be done by hand: `modules empty class argument`,
`modules TYPE_EMPTY_P`, `modules wrong-code parameter passing`,
`modules ABI empty`, plus the open wrong-code bugs that block PR103524.
