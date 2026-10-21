# Summary: [modules] ICE in finish_member_declaration instantiating V<int> via hidden friend
Component: c++   Keywords: ice-on-valid-code
Known to fail: 14.4.0, 15.2.0, 15.3.0, 16.2.0, 17.0 (trunk 20261020)   Known to work: (none found)
Blocks: c++-modules (PR103524)

## Description
Module `a` exports a class template `V` with a hidden-friend `operator*` and
an alias `V1 = V<int>`. Module `b` re-exports `a` and instantiates `V<int>` in
an inline function. An importer of `b` whose first completion of `V<int>`
happens inside a function body, through the friend operator (`2 * r` with
`r` a `V1&`), crashes in `finish_member_declaration`. Expected: the program
compiles and returns 0.

## Testcase
a.cppm
```c++
export module a;
export template <typename T> struct V {
  T x;
  friend V operator*(T s, V v) { return {s * v.x}; }
};
export using V1 = V<int>;
```

b.cppm
```c++
export module b;
export import a;
export inline V1 f(V1 v) { return v; }
```

main.cpp
```c++
import b;
void g(V1& r) { r = 2 * r; }
int main() { V1 v{1}; g(v); return v.x == 2 ? 0 : 1; }
```

## Command lines
```
g++ -std=c++26 -fmodules -Wall -Wextra -c a.cppm
g++ -std=c++26 -fmodules -Wall -Wextra -c b.cppm
g++ -std=c++26 -fmodules -Wall -Wextra -c main.cpp
```
<gcc -v output here>

## Compiler output
```
In module a, imported at b.cppm:2,
of module b, imported at main.cpp:1:
a.cppm: In instantiation of 'struct V@a<int>':
main.cpp:2:25:   required from here
    2 | void g(V1& r) { r = 2 * r; }
      |                         ^
a.cppm:2:39: internal compiler error: in finish_member_declaration, at cp/semantics.cc:4261
    2 | export template <typename T> struct V {
      |                                       ^
0x7f904e fancy_abort(char const*, int, char const*)
0xa0f259 instantiate_class_template(tree_node*)
0x82b8a9 build_new_op(op_location_t const&, tree_code, int, tree_node*, tree_node*, tree_node*, tree_node*, tree_node**, int)
0xa68f42 build_x_binary_op(op_location_t const&, tree_code, tree_node*, tree_code, tree_node*, tree_code, tree_node*, tree_node**, int)
0x9c847d c_parse_file()
0xaf5e79 c_common_parse_file()
```
The failing assertion is
`gcc_assert (TYPE_BEING_DEFINED (current_class_type) || LAMBDA_TYPE_P (TREE_TYPE (decl)))`.
Trunk 20261020 fails the same assertion at cp/semantics.cc:4324; 15.2 at
semantics.cc:4210; 14.4 at semantics.cc:3941.

## Notes
- No ICE if `V<int>` is spelled directly instead of through the exported
  alias `V1`.
- No ICE if `V1` is completed at namespace scope before the function, e.g.
  `static_assert(sizeof(V1) == sizeof(int));` after the import (the
  workaround in use).
- The same ICE happens when the using unit is an implementation unit of `b`
  (`module b;`) instead of a plain importer.
- Original code: `Vec3` = `BasicVec3<double>` defined in one module,
  instantiated in a second through a hidden friend in inline functions; the
  ICE came from `b[i] -= factor * b[k]` on a `std::vector<Vec3>` in an
  implementation unit of a third.
- Probably related: a "recursive lazy load" / "failed to load pendings"
  rejects-valid bug with the same ingredients (exported alias to a
  specialization that a second module instantiates), filed separately.
  <add PR number once filed>

## Duplicate check   (remove before filing)
Web search only (Bugzilla is behind a bot check):
`gcc bugzilla modules ICE "finish_member_declaration" hidden friend`,
`gcc "[modules]" ICE "finish_member_declaration" instantiate_class_template`,
`gcc bugzilla modules alias template specialization instantiated in another module hidden friend ICE`.
Only old non-modules ICEs in finish_member_declaration (PR16706, PR48029,
PR60390, PR67021, PR81917, PR84314) and PR123401 (friend `operator==<>`, not
modules). No duplicate found.

Final Bugzilla search to be done by hand: `finish_member_declaration modules`,
`[modules] ICE instantiate_class_template hidden friend`,
`modules alias specialization ICE`.
