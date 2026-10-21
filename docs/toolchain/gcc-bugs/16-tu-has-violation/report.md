# Summary: [modules] undefined reference to __tu_has_violation from imported inline function with pre()
Component: c++   Keywords: link-failure, rejects-valid
Known to fail: 16.1.0, 16.2.0, 17.0 (trunk 20261020)   Known to work: (none known)
Blocks: c++-modules (PR103524)

## Description
A module exports a non-template inline function with a precondition. A
translation unit that imports the module and calls the function, but contains
no contract of its own, fails to link. The inline function is emitted in the
importer, and its check calls `__tu_has_violation(...)`. GCC emits that helper
only in a translation unit that parses a contract itself, so the importer never
defines it. Expected: the program links and runs (exit code 0).

Adding any function with a contract to main.cpp makes the -O0 build link, but
not an -O1 or -O2 build (see Notes).

## Testcase
m.cppm:
```c++
export module m;
export inline int f(int x) pre(x > 0) { return x; }
```

main.cpp:
```c++
import m;
int main(int argc, char **) { return f(argc) - 1; }
```

## Command lines
```
g++ -std=c++26 -fmodules -Wall -Wextra -c m.cppm -o m.o
g++ -std=c++26 -fmodules -Wall -Wextra -c main.cpp -o main.o
g++ m.o main.o -lstdc++exp
```
<gcc -v output here>

## Compiler output
```
/usr/bin/ld: main.o: in function `f@m(int)':
m.cppm:2:(.text._ZW1m1fi[_ZW1m1fi]+0x1c): undefined reference to `__tu_has_violation(__builtin_contract_violation_type const&, unsigned short)'
collect2: error: ld returned 1 exit status
```
At -O2 the reference comes from `.text.startup`, that is from `main` after inlining.

## Notes
- Only non-template inline functions are affected. A function template with a
  precondition links, but its precondition is not checked in the importer
  (reported separately).
- If main.cpp parses a contract of its own, for example
  `[[gnu::used]] static void anchor(bool b) pre(b) {}`, the program links at
  -O0. It still fails at -O1 and -O2 with the same undefined reference: the
  importer's own helper is then inlined or cloned (`.constprop`), so no
  out-of-line `__tu_has_violation` exists for the imported function to call.
- `-fcontract-disable-optimized-checks -fcontract-checks-outlined` does not fix
  the testcase as written. Together with the anchor function above, the program
  links, but every contract-bearing function in the module unit then warns:
  ```
  m.cppm:2:19: warning: 'optimize' attribute currently unsupported in modules
      2 | export inline int f(int x) pre(x > 0) { return x; }
        |                   ^
  ```
  This warning has no `-W` option (it was added for PR108080), so `-Werror`
  builds cannot use these options. The user never wrote an `optimize`
  attribute; the contract options add one implicitly.
- Workaround: compile one contract in every importing unit and build importers
  at -O0, or keep contracts off inline functions in module interfaces.

## Duplicate check   (remove before filing)
- Bugzilla is behind an Anubis check, so web search only:
  "__tu_has_violation" (no hits), "gcc bugzilla contracts modules undefined
  reference", "contract_violation modules 2026", "'optimize' attribute
  currently unsupported in modules".
- Related, not the same: PR108205 ("Improve contracts support in modules", for
  the old contracts TS implementation); PR108080 and PR121396 (origin of the
  'optimize' warning).
- No duplicate found.
- Final Bugzilla search to be done by hand: `__tu_has_violation`;
  `contract undefined reference module`; component c++ bugs mentioning
  "contract" from 2026; `optimize attribute currently unsupported in modules`.
