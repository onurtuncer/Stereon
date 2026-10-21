# Summary: [modules] conflicting declaration 'struct std::contracts::contract_violation' with <source_location> in GMF
Component: c++   Keywords: rejects-valid, ice-on-valid-code
Known to fail: 16.1.0, 16.2.0 (error), 17.0 (trunk 20261020: ICE in mangle_decl)   Known to work: (none known)
Blocks: c++-modules (PR103524)

## Description
A module includes `<source_location>` in its global module fragment. That
header declares `namespace std::contracts { class contract_violation; }`. An
importing translation unit that contains a contract (so the compiler creates
its built-in `std::contracts::contract_violation`) and instantiates a template
from the module that uses `std::source_location` is rejected with a conflict
between the built-in type and the imported declaration. Expected: the two
should be merged (or the built-in reused), and the program should compile.
On trunk the diagnostic is gone, but the compiler segfaults in `mangle_decl`.

## Testcase
m.cppm:
```c++
module;
#include <source_location>
export module m;
export template <class T> int f() { return std::source_location::current().line(); }
```

main.cpp:
```c++
import m;
void h(int x) pre(x > 0) {}
int main() { return f<int>() - 4; }
```

## Command lines
```
g++ -std=c++26 -fmodules -Wall -Wextra -c m.cppm -o m.o
g++ -std=c++26 -fmodules -Wall -Wextra -c main.cpp -o main.o
g++ m.o main.o -lstdc++exp
```
<gcc -v output here>

## Compiler output
16.2.0, compiling main.cpp:
```
<built-in>: error: conflicting declaration 'struct std::contracts::contract_violation'
In file included from m.cppm:2,
of module m, imported at main.cpp:1:
/usr/include/c++/16.2.0/source_location:42:31: note: previous declaration as 'class std::contracts::contract_violation'
   42 |   namespace contracts { class contract_violation; }
      |                               ^~~~~~~~~~~~~~~~~~
```
(Replace the header path with your local one.)

trunk 20261020, compiling main.cpp:
```
main.cpp:3:35: internal compiler error: Segmentation fault
    3 | int main() { return f<int>() - 4; }
      |                                   ^
0x2c5c7fb internal_error(char const*, ...)
0xd0f079 mangle_decl(tree_node*)
0x1996945 decl_assembler_name(tree_node*)
0x104973e symbol_table::insert_to_assembler_name_hash(symtab_node*, bool)
0x1067936 symbol_table::finalize_compilation_unit()
```

## Notes
- Without the contract in main.cpp (remove `h`), the program compiles and
  runs on 16.2 and on trunk.
- No library code is needed. The same error (and the same trunk ICE) occurs
  with the header replaced by a plain declaration in the global module
  fragment. That version warns `-Wglobal-module` and declares a name in `std`,
  so the testcase above uses the real header:
  ```c++
  module;
  namespace std::contracts { class contract_violation; }
  export module m;
  export template <class T> int f() { std::contracts::contract_violation *p = nullptr; return p != nullptr; }
  ```
- It was first seen with `<meta>` in the global module fragment and a
  reflection-based serialiser instantiated in the importer (`-freflection`);
  `<meta>` includes `<source_location>`. Reflection is not needed.
- The class-key difference (`struct` built-in vs `class` in libstdc++) is not
  the issue in itself. The built-in declaration is not merged with the
  imported one, so this belongs to the c++ component, not libstdc++.
- Workaround: instantiate such templates only in units that contain no contract.

## Duplicate check   (remove before filing)
- Web search only (Bugzilla is behind Anubis): "contract_violation conflicting
  declaration modules", "gcc bugzilla contracts modules 2026".
- Related, not the same: PR108205 (old contracts TS: generated
  contract_violation type and built-in declarations vs modules; same class of
  problem for the P2900 implementation); PR125787 ([modules] conflicting
  imported declarations with -freflection; no contracts involved).
- Final Bugzilla search to be done by hand: `contract_violation`;
  `conflicting declaration contracts`; `mangle_decl module contract`.
