# Summary: [modules] pre/post of function templates imported from a module are not checked
Component: c++   Keywords: wrong-code
Known to fail: 16.1.0, 16.2.0, 17.0 (trunk 20261020)   Known to work: (none known)
Blocks: c++-modules (PR103524)

## Description
A precondition (or postcondition) on a function template exported from a
module is silently not evaluated when the template is instantiated and called
in an importing translation unit. The program below violates the precondition
and should report the violation and terminate (enforce is the default
semantic). Instead it exits normally with status 0. The same template
without modules reports the violation.

## Testcase
m.cppm:
```c++
export module m;
export template <class T> T g(T x) pre(x > 0) { return x; }
```

main.cpp:
```c++
import m;
int main(int argc, char **) { return g(argc - 1); }
```

## Command lines
```
g++ -std=c++26 -fmodules -Wall -Wextra -c m.cppm -o m.o
g++ -std=c++26 -fmodules -Wall -Wextra -c main.cpp -o main.o
g++ m.o main.o -lstdc++exp
./a.out; echo $?
```
<gcc -v output here>

## Compiler output
No diagnostics. Running prints nothing and exits with status 0. (The program
even links without `-lstdc++exp`, so no check is emitted at all.)

Expected, as the non-module version prints (the template defined in main.cpp
instead of being imported):
```
contract violation in function T g(T) [with T = int] at main.cpp:1: x > 0
[assertion_kind: pre, semantic: enforce, mode: predicate_false, terminating: yes]
terminate called without an active exception
```
(SIGABRT, exit status 134)

## Notes
- Same result at -O0 and -O2, with `-fcontracts-client-check=pre`, and when
  main.cpp contains contracts of its own.
- `post(r : r > 0)` on an imported template is ignored the same way, and so is
  `pre` on a member function of an imported class template
  (`template <class T> struct A { T get(T x) const pre(x > 0); };`).
- A `contract_assert` inside the body of an imported function template is
  checked: `export template <class T> T a(T x) { contract_assert(x > 0); return x; }`
  reports the violation in the importer. Only the function contract
  specifiers seem to be lost.
- Contracts on non-template inline functions from the module are checked in
  the importer, but they hit a link error unless the importer has a contract
  of its own (reported separately: undefined reference to `__tu_has_violation`).
- No workaround: the checks fire only for calls instantiated inside the
  module's own units.

## Duplicate check   (remove before filing)
- Web search only (Bugzilla is behind Anubis): "gcc contracts modules
  precondition template imported not checked", "gcc bugzilla 2026 contracts
  modules". No matching report. Related only: PR108205 (old contracts TS + modules).
- Final Bugzilla search to be done by hand: `contract template module`;
  `pre module import ignored`; component c++ bugs mentioning "contract" from
  2026; the PR103524 dependency list.
