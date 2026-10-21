# NOT TO BE FILED: fixed on trunk as of 20261020 (not a modules bug)
# Summary: "contract condition is not constant" for pre() calling a member already called in the same constant evaluation
Component: c++   Keywords: rejects-valid
Known to fail: 16.1.0, 16.2.0   Known to work: 17.0 (trunk 20261020)

## Description
A `static_assert` evaluates a constexpr function that calls `s.size()` and then
`s.get()` on a `const S&` parameter. `get()` has the precondition `size() > 0`,
which holds. GCC 16.2 rejects the constant evaluation with "contract condition
is not constant". If the function does not call `size()` before `get()`, or
works on a local object instead of a reference parameter, it compiles. Modules
are not involved. Trunk accepts the testcase and still diagnoses a false
predicate in the same shape ("contract predicate is false in constant
expression"), so the bug appears to be fixed there.

## Testcase
main.cpp:
```c++
struct S {
  int n = 1;
  constexpr int size() const { return n; }
  constexpr int get() const pre(size() > 0) { return n; }
};

constexpr int f(const S &s) { return s.size() + s.get(); }

static_assert(f(S{}) == 2);

int main() { return f(S{}) - 2; }
```

## Command lines
```
g++ -std=c++26 -Wall -Wextra main.cpp -lstdc++exp
```
<gcc -v output here>

## Compiler output
```
main.cpp:4:29: error: contract condition is not constant
    4 |   constexpr int get() const pre(size() > 0) { return n; }
      |                             ^~~~~~~~~~~~~~~
```

## Notes
- Reduced from Stereon's expansion arithmetic: `Expansion<N>`'s converting
  constructor did `: m_Length(e.GetLength())` and then read `e[i]`, whose
  precondition is `i < GetLength()`. The full code fails the same way as a
  single non-module file. With modules, the error appears in importers only
  because that is where the `static_assert`s were.
- Variations that compile on 16.2: no earlier `s.size()` call; `S` a local
  variable of the constexpr function; the predicate reading `n` directly
  instead of calling `size()`.
- Related observation on 16.2 (also fixed on trunk): preconditions of
  non-member constexpr functions are not evaluated during constant evaluation
  at all. `constexpr int g(int x) pre(x > 0) { return x; } static_assert(g(0) == 0);`
  is accepted by 16.2 and rejected by trunk.
- Rather than a new PR, consider asking for a backport to the gcc-16 branch
  once the fixing commit is identified.

## Duplicate check   (remove before filing)
- Web search for "contract condition is not constant" found no report. The
  message comes from `check_for_failed_contracts` in cp/constexpr.cc.
- Final Bugzilla search to be done by hand: `contract condition is not
  constant`; component c++ bugs mentioning "contract" and "constexpr" from
  2026 (to find the fixing PR for a backport request).
