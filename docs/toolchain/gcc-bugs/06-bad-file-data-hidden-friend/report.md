# Summary: [modules] "Bad file data" loading inline fn calling hidden friend across three modules
Component: c++   Keywords: rejects-valid
Known to fail: 13.4.0 (-fmodules-ts), 14.3.0 (-fmodules-ts), 15.2.0, 16.1.0, 16.2.0, 17.0 (trunk 20261021)   Known to work: none found
Blocks: c++-modules (PR103524)
See also: PR99389, PR99426, PR125635

## Description
Three modules form a chain a <- b <- c. Module `a` defines class templates
`V<T>` and `P<T>`. `P<T>` has a hidden friend `operator-` that returns `V<T>`.
`P<int>` is implicitly instantiated in `a`, and `V<int>` is implicitly
instantiated in `b`. Module `c` has an exported inline non-template function
that calls `P<int>`'s hidden friend `operator-`. All four units compile
without diagnostics until the importer: a unit that imports `c` and names `f`
cannot read `c.gcm` ("failed to read compiled module cluster N: Bad file
data", then "failed to load binding '::f@c'"). Expected: the program compiles
and exits with 0.

## Testcase
a.cppm
```c++
export module a;
export template <typename T> struct V { T x; };
export template <typename T> struct P {
  T x;
  friend V<T> operator-(const P& p, const P& q) { return {p.x - q.x}; }
};
export inline int g(const P<int>& p) { return p.x; }
```

b.cppm
```c++
export module b;
export import a;
export struct S { V<int> v; };
```

c.cppm
```c++
export module c;
export import b;
export inline int f(const P<int>& p, const P<int>& q) { return (p - q).x; }
```

main.cpp
```c++
import c;
int main() { P<int> p{1}; return f(p, p); }
```

## Command lines
```
g++ -std=c++26 -fmodules -Wall -Wextra -c a.cppm
g++ -std=c++26 -fmodules -Wall -Wextra -c b.cppm
g++ -std=c++26 -fmodules -Wall -Wextra -c c.cppm
g++ -std=c++26 -fmodules -Wall -Wextra -c main.cpp
```
<gcc -v output here>

## Compiler output
```
In module imported at main.cpp:1:1:
c: In function 'int main()':
c: error: failed to read compiled module cluster 1: Bad file data
c: note: compiled module file is 'gcm.cache/c.gcm'
c: error: failed to read compiled module cluster 2: Bad file data
main.cpp:2:34: fatal error: failed to load binding '::f@c'
    2 | int main() { P<int> p{1}; return f(p, p); }
      |                                  ^
compilation terminated.
```
Same at -O0 and -O2, with -std=c++20, and on a checking-enabled 16.2 build
(no assertion fires while `c.gcm` is written).

## Notes
Each of these changes makes the testcase build and run:
- remove `g` from `a`, so that `P<int>` is first instantiated in `c`;
- remove `S` from `b`, so that `V<int>` is not instantiated in `b`;
- make `operator-` a member function instead of a hidden friend;
- have the hidden friend return `T` instead of `V<T>`;
- merge `b` and `c` into one module.

If `main` itself evaluates `p - p` before calling `f`
(`return (p - p).x + f(p, p);`), the CMI is read without error, but the
assembler then reports `symbol '_ZW1amiRKS_1PIiES3_' is already defined`:
`P<int>`'s friend `operator-` is emitted twice in `main.o`. This may be the
same defect seen from the other side.

Found in a real project, where it showed up as "failed to load binding
'::sin@...'" in a unit that included `<cmath>` and imported the module; the
bad cluster happened to hold `::sin`'s binding, `<cmath>` is not involved.
There, `Point3 - Point3` is a hidden friend of a class template the first
module instantiates, the second module instantiates the vector type as a
member, and the third uses `Point3 - Point3` in an inline non-template
function.

Workaround: move such inline functions out of line (into an implementation
unit).

## Duplicate check   (remove before filing)
Web search and the gcc-bugs mail archive only (Bugzilla is behind a bot check):
"failed to read compiled module cluster" "Bad file data"; "failed to load
binding" "Bad file data"; "failed to load pendings" hidden friend.

Candidates:
- PR99389 "[modules] bad serialization of data" (RESOLVED FIXED, r11-7524):
  same symptom for an exported inline function using a member of a class
  template instantiation. Still fails on 16.2 and trunk, so a different path
  or an incomplete fix.
- PR99426 (open, same generic error text, header units / standard library).
- PR125635 (UNCONFIRMED, 2026-06): also a three-module chain, "failed to load
  pendings", std::list and inheritance. Possibly the same root cause.
- PR116805: same error text, different trigger.

No clear duplicate: file new, with See also 99389, 99426, 125635.
Final Bugzilla search to be done by hand: `"Bad file data" modules`,
`"failed to load binding"`, `"failed to read compiled module cluster"`,
`"failed to load pendings"`, `hidden friend modules instantiation`, and the
open dependencies of PR103524.
