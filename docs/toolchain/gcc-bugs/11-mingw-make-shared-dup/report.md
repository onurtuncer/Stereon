Merged into ../10-mingw-format-dup/report.md (testcase 2). Same mechanism: a
function-local static of a libstdc++ inline function, instantiated in an
importer through a module template, is emitted WEAK instead of UNIQUE on ELF
and gives a multiple definition on mingw. This directory holds testcase 2 and
its repro.sh.
