# Brief for preparing GCC bug reports

Goal: one upstream-ready GCC Bugzilla report per compiler bug listed in
`docs/toolchain/gcc16-feature-probe.md` ("Known compiler bugs"). Each report
gets its own directory `docs/toolchain/gcc-bugs/NN-slug/` holding:

- the testcase files: module interfaces `*.cppm`, implementation units
  `*.cpp`, and `main.cpp` as the LAST translation unit;
- `repro.sh`: the exact commands for a local GCC 16.2 (one `g++` call per
  TU, in order, then link/run if needed), so the user can re-check on MSYS2
  and Linux;
- `report.md`: the Bugzilla text (template below).

## Verifying on upstream GCC (no GCC 16 on this machine)

Use Compiler Explorer through `docs/toolchain/gcc-bugs/ce.ps1` (PowerShell):

    & docs\toolchain\gcc-bugs\ce.ps1 -Dir <dir> -Order a.cppm,b.cpp -Execute `
        -Compiler g162 -Flags "-std=c++26 -fmodules -fdiagnostics-color=never"

All non-main files in `<dir>` are uploaded; `-Order` lists those compiled
before `main.cpp`, in order (GCC writes each module's CMI to gcm.cache before
compiling the next file). Use `-Execute` to build and run,
`-CompileOnly` to only compile (no link). Compilers: `g162` (16.2), `g161`,
`gsnapshot` (trunk), `g162assert` (checking build: useful for ICE messages).
Add `-O2`, `-freflection`, `-fcontract-evaluation-semantic=enforce`,
`-lstdc++exp` etc. to `-Flags` as needed. CE runs Linux x86-64 glibc.
Be gentle with the API: no tight loops, a few seconds between calls is fine.

Only send small synthetic testcases or Stereon source (Stereon is public,
MPL-2.0). Never send anything else.

## Reduction rules (from https://gcc.gnu.org/bugs/)

- Minimal: no Stereon names left, no unnecessary includes; ideally no
  `#include` at all. Fewest files and lines that still show the exact bug.
- Must reproduce exactly: no `...` placeholders in testcases.
- Compile cleanly with `-Wall -Wextra` apart from the bug itself.
- If a bug does NOT reproduce on upstream 16.2 Linux (e.g. MinGW-only), say so
  plainly in the report and in your answer; don't fake a reproduction. A
  MinGW-only bug may still be reported but must say which target it was seen
  on and that it was not reproduced on Linux.
- Record whether trunk (`gsnapshot`) still fails. If trunk is fixed, say so:
  that bug is NOT to be filed (note "fixed on trunk as of <date>").

## Duplicate check

Search GCC Bugzilla before writing, e.g. with WebFetch on
`https://gcc.gnu.org/bugzilla/buglist.cgi?quicksearch=<terms>&component=c%2B%2B`
(also try `component=libstdc%2B%2B`), using the key error text
("Bad file data", "recursive lazy load", "mangle_module", "TU-local", ...)
and the `modules` keyword/`c++-modules` meta bug. Open plausible candidates.
Record what you searched and any likely duplicates in report.md (a
"Duplicate check" section, which the user removes before filing). If a clear
duplicate exists, write the report as a comment to add to that PR instead.

## report.md template

    # Summary: [modules] <one line, as Bugzilla Summary; <= 90 chars>
    Component: c++ (or libstdc++)   Keywords: ice-on-valid-code / rejects-valid / wrong-code / link-failure
    Known to fail: 16.1.0, 16.2.0, 17.0 (trunk YYYYMMDD)   Known to work: (if any)
    Blocks: c++-modules (PR103524)

    ## Description
    2-6 sentences: what happens, what was expected.

    ## Testcase
    Each file in a fenced block, with its name.

    ## Command lines
    The exact g++ commands, one per TU.

    ## Compiler output
    Verbatim (trim backtrace to the first ~10 frames).

    ## Notes
    Workaround, observations (which variations make it go away).

    ## Duplicate check   (remove before filing)
    Searches and candidates.

The user will add `gcc -v` output themselves; leave a line
`<gcc -v output here>` under "Command lines".
