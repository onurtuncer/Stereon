# GCC 16.2 bug reports (modules, contracts)

Prepared 2026-10-21 from the "Known compiler bugs" table in
`../gcc16-feature-probe.md`. Each directory holds a reduced testcase, a
`repro.sh` (one `g++` call per unit) and `report.md` (Bugzilla text).
Reproductions were checked on upstream x86_64-linux-gnu GCC through Compiler
Explorer (`ce.ps1`), and 10 on its MinGW 16.2.

| # | Bug | Kind | Fails | Action |
| --- | --- | --- | --- | --- |
| 01 | defaulted hidden-friend `operator==` from a module | ICE | 16.1, 16.2, trunk; OK on 14.1, 15.2 | **File** (16 regression) |
| 02 | stdexec under `-fmodules` (extern then inline variable) | ICE | 16.1, 16.2 | Don't file: PR126783, fixed for 16.3 |
| 03 | `std::source_location` "TU-local" | rejects-valid | not reproduced | Re-check locally first |
| 04 | "contract condition is not constant" | rejects-valid | 16.1, 16.2; fixed on trunk | Don't file; maybe ask for a backport |
| 05 | "Bad file data" (Cylinder / BasicCircle) | rejects-valid | — | Same as 06 |
| 06 | "Bad file data" via hidden friend across three modules | rejects-valid | 13.4 to trunk | **File** |
| 07 | empty-class argument before memory argument, x86-64 | **wrong-code** | 12.2 to trunk | **File** |
| 08 | "recursive lazy load" via alias to a specialization | rejects-valid | 15.3, 16.x, trunk; OK on 15.2 | **File** (regression) |
| 09 | inline explicit specialization of a variable template | link-failure | 14.4 to trunk | **File** |
| 10 | local static of libstdc++ inline function, mingw | link-failure | 16.2 mingw | **File** after local asm check (includes 11) |
| 11 | (`make_shared` variant of 10) | | | Merged into 10 |
| 12 | conflicting `std::contracts::contract_violation` | rejects-valid / ICE on trunk | 16.1, 16.2, trunk | **File** |
| 13 | "Bad file data" from inline partition functions | rejects-valid | — | Same as 06 |
| 15 | `finish_member_declaration` ICE via alias + hidden friend | ICE | 14.4 to trunk | **File**, cross-reference 08 |
| 16 | undefined `__tu_has_violation` | link-failure | 16.1, 16.2, trunk | **File** |
| 17 | `pre`/`post` of imported templates not checked | **wrong-code** | 16.1, 16.2, trunk | **File** |

(14 is unused: the probe table's third "Bad file data" row went into 13.)

## Before filing each report
1. Run `repro.sh` with the local GCC 16.2 (MSYS2 and, if possible, Linux) and
   paste `gcc -v` over the `<gcc -v output here>` line.
2. Search Bugzilla by hand with the terms in the report's "Duplicate check"
   section (Bugzilla blocks automated searches), then delete that section.
3. Attach the testcase files, or keep them inline as in the report.
4. For 08/15: file one, then put its PR number in the other.

## Corrections the probe table needs
- 07: the cause is the empty tag argument (`Key{}`) before a by-value
  aggregate, not the `string_view` aggregate. Other module functions taking
  an empty class by value before a large by-value argument are miscompiled
  on Linux too.
- 06/05/13: the trigger is a hidden friend of a class template used in an
  inline function of a third module; `<cmath>` and `-O2` are not involved.
- 03: libstdc++ does not declare `source_location` with internal linkage.
- 12: `<source_location>` is enough; `<meta>` and `-freflection` are not needed.
- 17: `contract_assert` in an imported template body is checked; only
  `pre`/`post` are dropped.
- 02: fixed in 16.3 (PR126783); `stereon_exec_sources()` can go then.
