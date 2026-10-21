# Summary: [15/16/17 Regression] [modules] "recursive lazy load" importing a partition that instantiates A<V1>
Component: c++   Keywords: rejects-valid
Known to fail: 15.3.0, 16.1.0, 16.2.0, 17.0 (trunk 20261020)   Known to work: 14.4.0, 15.1.0, 15.2.0
Blocks: c++-modules (PR103524)

## Description
Module `a` exports a class template `V`, an alias `V1 = V<int>` and a
second class template `A`. Module `b` re-exports `a` and instantiates `V<int>`
in an inline function. A partition `c:part` imports `b` and instantiates
`A<V1>` in a function template. Compiling the primary interface `c`, which only
does `export import :part;`, fails with "recursive lazy load" / "failed to load
pendings for '::V'". The code is valid and compiles with GCC 14.4, 15.1 and
15.2.

## Testcase
a.cppm
```c++
export module a;
export template <typename T> struct V { T x; };
export using V1 = V<int>;
export template <typename T> struct A { T e[4]; };
```

b.cppm
```c++
export module b;
export import a;
export inline V1 f(V1 v) { return v; }
```

c-part.cppm
```c++
export module c:part;
import b;
export template <typename S>
V1 g(S) {
  A<V1> a{};
  return a.e[2];
}
```

c.cppm
```c++
export module c;
export import :part;
```

main.cpp
```c++
import c;
int main() {}
```

## Command lines
```
g++ -std=c++26 -fmodules -Wall -Wextra -c a.cppm
g++ -std=c++26 -fmodules -Wall -Wextra -c b.cppm
g++ -std=c++26 -fmodules -Wall -Wextra -c c-part.cppm
g++ -std=c++26 -fmodules -Wall -Wextra -c c.cppm
g++ -std=c++26 -fmodules -Wall -Wextra -c main.cpp
```
<gcc -v output here>

## Compiler output
```
c.cppm:2:21: error: recursive lazy load
    2 | export import :part;
      |                     ^
c.cppm:2:21: fatal error: failed to load pendings for '::V'
compilation terminated.
```
(main.cpp then fails because gcm.cache/c.gcm was not written.) The
checking-enabled 16.2 build gives the same diagnostic, with no ICE.

## Notes
- Compiles with `-fno-module-lazy`.
- Compiles if `V<int>` is spelled directly everywhere instead of through the
  exported alias `V1`.
- Compiles if module `b` is removed and `c:part` imports `a` directly, so the
  failure needs `V<int>` to be instantiated in a module other than the one
  that defines `V`.
- Compiles if `c:part` is an ordinary module imported by main.cpp. The failure
  happens when the partition is imported into its primary interface.
- `g` does not need to be a template: an inline non-template function fails
  the same way.
- Original code: a function template in a partition held a local
  `std::array<Vec3, 4>`, where `Vec3` is an alias of `BasicVec3<double>` from
  another module; every implementation unit of the module then failed with
  "failed to load pendings for 'Stereon::BasicVec3'". Workaround: name the
  element type through a dependent `decltype`.
- On Compiler Explorer it works with 15.2 and fails with 15.3, 16.1, 16.2 and
  trunk, so the cause is probably a change backported to the 15 branch.
- Probably related: an ICE in finish_member_declaration with the same
  ingredients (exported alias to a specialization that a second module
  instantiates), filed separately. <add PR number once filed>

## Duplicate check   (remove before filing)
Web search only (Bugzilla is behind a bot check):
`gcc modules "recursive lazy load" "failed to load pendings"`,
`"recursive lazy load" gcc modules bug` (gcc.gnu.org).

Candidate: PR125635 "[modules] failed to load pendings for
std::__cxx11::list after unused module import" (libstdc++, 16.0, UNCONFIRMED).
Similar shape (a type instantiated across three modules), but its first error
is "Bad file data", not "recursive lazy load", and its testcase needs
inheritance from a polymorphic class. Possibly the same root cause;
cross-reference it.

Related history: PR99377, PR116317, PR99426.

Final Bugzilla search to be done by hand: `"recursive lazy load"`,
`"failed to load pendings" modules`, `modules alias specialization lazy`.
