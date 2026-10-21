# Summary: [modules] local static of inline function instantiated in importer: multiple definition on mingw
Component: c++   Keywords: link-failure
Known to fail: 16.2.0 (x86_64-w64-mingw32)   Known to work: (none known)
Blocks: c++-modules (PR103524)

## Description
A module exports a function template whose body calls a libstdc++ inline
function with a function-local static:
`std::__format::__write_escaped_unicode_part<...>::__replace_rep` (from
`std::format` with a `string_view` argument) or
`std::_Sp_make_shared_tag::_S_ti()::__tag` (from `std::make_shared`). An
importer instantiates the template. Another TU uses the same libstdc++
function through an ordinary `#include`. On x86_64-w64-mingw32 the link fails
with "multiple definition" of the static. The two copies should fold like any
other vague-linkage local static.

On x86_64-linux-gnu the program links, but `readelf -sW` shows the importer's
copy bound `WEAK`, while the copy from the `#include` TU is `UNIQUE`
(gnu_unique_object). The importer path therefore emits the static with
different linkage. On PE-COFF that presumably means no `.linkonce` section, so
the importer's copy becomes a hard duplicate.

## Testcase 1 (std::format)
m.cppm:
```c++
module;
#include <format>
#include <string_view>
export module m;
export template <class T> unsigned fmt(const T* s) { return std::format("{}", std::string_view(s)).size(); }
```
a.cpp:
```c++
import m;
unsigned fa() { return fmt("x"); }
```
main.cpp:
```c++
#include <format>
#include <string_view>
unsigned fa();
int main() { return std::format("{}", std::string_view("y")).size() + fa() - 2; }
```

## Testcase 2 (std::make_shared)
m.cppm:
```c++
module;
#include <memory>
export module m;
export template <class T> std::shared_ptr<T> make(T v) { return std::make_shared<T>(v); }
```
a.cpp:
```c++
#include <typeinfo>
import m;
int fa() { return *make(1); }
```
main.cpp:
```c++
#include <memory>
int fa();
int main() { return *std::make_shared<int>(1) + fa() - 2; }
```
(Without `#include <typeinfo>` in a.cpp, g++ rejects the instantiation of
`_Sp_counted_ptr_inplace::_M_get_deleter` in the importer with "must
'#include <typeinfo>' before using 'typeid'". That may be worth a separate
look.)

## Command lines (both testcases)
```
g++ -std=c++26 -fmodules -x c++ -c m.cppm -o m.o
g++ -std=c++26 -fmodules -c a.cpp -o a.o
g++ -std=c++26 -fmodules -c main.cpp -o main.o
g++ m.o a.o main.o
```
<gcc -v output here>

## Compiler output (MinGW-w64 GCC 16.2.0, posix, ucrt)
(Paths shortened from the Compiler Explorer run; replace with the verbatim
MSYS2 output before filing.)
Testcase 1:
```
ld.exe: main.o:.../include/c++/16.2.0/format:1155: multiple definition of `std::__format::__write_escaped_unicode_part<char, std::__format::_Sink_iter<char> >(std::__format::_Sink_iter<char>, std::basic_string_view<char, std::char_traits<char> >&, bool&, std::__format::_Term_char)::__replace_rep'; a.o:.../format:1155: first defined here
collect2.exe: error: ld returned 1 exit status
```
Testcase 2:
```
ld.exe: main.o:.../bits/shared_ptr_base.h:633: multiple definition of `std::_Sp_make_shared_tag::_S_ti()::__tag'; a.o:.../shared_ptr_base.h:633: first defined here
collect2.exe: error: ld returned 1 exit status
```

## Notes
- x86_64-linux-gnu GCC 16.2: links. With main.cpp reduced to
  `int fa(); int main(){return fa()-1;}`, `readelf -sW a.o` shows the static
  as `OBJECT WEAK` in both testcases. If a.cpp #includes `<format>` /
  `<memory>` itself and calls the function directly, it is `OBJECT UNIQUE`.
- On mingw the #include TU emits
  `.section .rdata$_ZZNSt19_Sp_make_shared_tag5_S_tiEvE5__tag,"dr"` with
  `.linkonce same_size`.
- If the static is emitted by a non-template function in the module interface
  (`export int fa() { return *std::make_shared<int>(1); }`), it also gets
  `.linkonce same_size` and the program links. The failure needs the
  instantiation to happen in the importer.
- Workaround: the importer also includes the libstdc++ header in its global
  module fragment; the copies then fold.
- Not checked on trunk (no mingw trunk compiler available).

## Duplicate check   (remove before filing)
- Before filing, run repro.sh on MSYS2 and confirm the a.s / main.s
  difference (`.linkonce` missing in a.s). Compiler Explorer cannot emit asm
  for an importer TU, so this is not yet confirmed.
- Web search found only PR115686 (mingw inline thread_local TLS init, no
  modules) and PR90562, both unrelated. A gcc-mirror commit search for
  "modules gnu_unique" / "modules local static comdat" found nothing.
- Final Bugzilla search to be done by hand: `modules multiple definition mingw`,
  `modules gnu_unique`, `modules local static weak`, `__replace_rep`,
  `_Sp_make_shared_tag`.
