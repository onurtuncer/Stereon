# Summary: [modules] undefined reference to inline explicit specialization of variable template
Component: c++   Keywords: link-failure
Known to fail: 14.4.0, 15.3.0, 16.2.0, 17.0 (trunk 20261020)   Known to work: (none found)
Blocks: c++-modules (PR103524)

## Description
A module interface exports a variable template `N` and declares an `inline
constexpr` explicit specialization `N<L>` (outside any export block). An
importer that odr-uses `N<L>`, here by taking its address, references the
symbol `N@a<L@a>`. Neither object file defines it, so linking fails. Expected:
the program links and returns 0; the inline variable should be emitted
(COMDAT) wherever it is odr-used, or by the module's object.

## Testcase
a.cppm
```c++
export module a;
export struct L {};
export template <typename T> inline constexpr int N = 0;
template <> inline constexpr int N<L> = 1;
```

main.cpp
```c++
import a;
int main() {
  const int* p = &N<L>;
  return *p - 1;
}
```

## Command lines
```
g++ -std=c++26 -fmodules -Wall -Wextra -c a.cppm
g++ -std=c++26 -fmodules -Wall -Wextra -c main.cpp
g++ a.o main.o -o main
```
<gcc -v output here>

## Compiler output
```
ld: main.o: in function `main':
main.cpp:3:(.text+0x8): undefined reference to `N@a<L@a>'
collect2: error: ld returned 1 exit status
```

## Notes
- Links if the specialization is not `inline`
  (`template <> constexpr int N<L> = 1;`): a.o then defines it.
- Links if the interface itself odr-uses the specialization, for example
  `export const int* use() { return &N<L>; }` in a.cppm.
- Reading `N<L>` only in constant expressions avoids the odr-use and works.
- Seen originally with
  `template <> inline constexpr std::string_view GeometryTypeName<Line> = "Line";`
  in a module partition, odr-used by `std::format` in another module and in
  tests.

## Duplicate check   (remove before filing)
Web search only (Bugzilla is behind a bot check):
`gcc modules explicit specialization variable template "undefined reference" inline`
(gcc.gnu.org), `"[modules]" variable template specialization not emitted importer undefined reference gcc bug`.
Only non-modules variable-template link bugs (PR70056, PR70095, PR69515,
PR88417). No duplicate found.

Final Bugzilla search to be done by hand: `modules variable template
specialization undefined reference`, `modules inline variable not emitted`,
`[modules] explicit specialization inline variable`.
