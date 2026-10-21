# Summary: [16/17 Regression] [modules] ICE in mangle_module with imported defaulted friend operator==
Component: c++   Keywords: ice-on-valid-code
Known to fail: 16.1.0, 16.2.0, 17.0 (trunk 20261021)   Known to work: 14.1.0, 15.2.0
Blocks: c++-modules (PR103524)

## Description
A module exports a class with a defaulted hidden-friend `operator==`. A
translation unit that imports the module and compares two objects with `==`
crashes the compiler while mangling the implicitly defined operator. GCC 15.2
compiles and runs the same code correctly.

The ICE needs the comparison in the importer: the module unit alone compiles,
and an importer that does not use `==` compiles. A defaulted member
`bool operator==(const S&) const = default;` works. The crash is the same with
`-std=c++20` and `-std=c++26`, with and without `constexpr`.

## Testcase
m.cppm:
```c++
export module m;
export struct S {
  int x;
  friend bool operator==(const S&, const S&) = default;
};
```

main.cpp:
```c++
import m;
int main() { S a{1}, b{1}; return a == b ? 0 : 1; }
```

## Command lines
```
g++ -std=c++20 -fmodules -c m.cppm
g++ -std=c++20 -fmodules -c main.cpp
```
<gcc -v output here>

## Compiler output
Release build of 16.2.0 (x86_64-linux-gnu):
```
In module m, imported at main.cpp:1:
m.cppm:4:15: internal compiler error: Segmentation fault
    4 |   friend bool operator==(const S&, const S&) = default;
      |               ^~~~~~~~
0x91c395 module_state::mangle(bool)
0x9004c1 mangle_decl(tree_node*)
0x13a2401 decl_assembler_name(tree_node*)
0xbb0bc6 symbol_table::finalize_compilation_unit()
```
Checking build of 16.2.0:
```
m.cppm:4:15: internal compiler error: in mangle_module, at cp/module.cc:16786
```
Trunk (20261021):
```
m.cppm:4:25: internal compiler error: in mangle_module, at cp/module.cc:16810
```

## Notes
Workaround: declare the defaulted comparison as a member.

Seen first on x86_64-w64-mingw32 (MSYS2 UCRT64 GCC 16.2.0) and also with
partitions, `-freflection` and contracts. Reproduced above on upstream
x86_64-linux-gnu builds.

## Duplicate check   (remove before filing)
- Web search for "mangle_module" / "module_state::mangle" + defaulted friend
  operator==: no GCC report found (closest: PR106361, PR109790, PR123401,
  none modules-related).
- Final Bugzilla search to be done by hand: `mangle_module`,
  `modules defaulted friend`, `modules operator== default ICE`, and the
  dependency list of PR103524.
