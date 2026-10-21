# DO NOT FILE: duplicate of PR c++/126783, fixed on trunk and releases/gcc-16

Fixed on trunk 2026-08-13 (a03b3ab4b0a6, "c++/modules: ICE with GMF variable
later defined inline [PR126783]") and backported to releases/gcc-16
(61f16b5960e4), so it should be fixed in 16.3. Testcase g++.dg/modules/gmf-6.C
has the same pattern. Kept for Stereon's notes: recheck with 16.3, then drop
the `stereon_exec_sources()` workaround.

Verified 2026-10-21 on Compiler Explorer: the testcase ICEs on 16.1.0 and
16.2.0, builds on 15.2.0 and on trunk 20261020.

## Testcase
main.cpp:
```c++
extern const int i;
inline const int i = 0;
int main() { return i; }
```
Also fails with `extern int i; inline int i = 0;` and with a class type
(`inline constexpr S s{};`). No module declaration or import is needed.

## Command line
```
g++ -std=c++26 -fmodules -c main.cpp
```

## Compiler output (16.2.0)
```
main.cpp:2:18: internal compiler error: Segmentation fault
    2 | inline const int i = 0;
      |                  ^
0x93834e transfer_defining_module(tree_node*, tree_node*)
0x8996b6 duplicate_decls(tree_node*, tree_node*, bool, bool)
0x957e0c pushdecl(tree_node*, bool)
```
Checking build: `internal compiler error: in transfer_defining_module, at cp/module.cc:22398`.

## How stdexec hits it
stdexec declares its CPOs `extern` first (`extern forwarding_query_t const
forwarding_query;` in `__execution_fwd.hpp`) and defines them `inline
constexpr` later (`__query.hpp`). So every TU that includes stdexec and is
compiled with `-fmodules` ICEs, and CMake adds `-fmodules` to every source it
scans. Without `-fmodules` the file compiles cleanly.

Not investigated: the separate stdexec header-unit error ("no class template
named '__f' in 'struct stdexec::just_t'").

## Duplicate check
Found through a gcc-mirror commit search for `transfer_defining_module`.
By hand: open PR126783 and confirm its target milestone is 16.3.
