#!/bin/sh
# GCC 16.2: "conflicting declaration 'struct std::contracts::contract_violation'"
# when a unit that parses a contract instantiates a template from a module
# whose global module fragment includes <source_location> (or <meta>).
# Expected: compiles, links, exits 0. Trunk 20261020 ICEs in mangle_decl instead.
set -x
CXX=${CXX:-g++}
FLAGS="-std=c++26 -fmodules -Wall -Wextra"
rm -rf gcm.cache *.o a.out

$CXX $FLAGS -c m.cppm -o m.o
$CXX $FLAGS -c main.cpp -o main.o
$CXX m.o main.o -lstdc++exp -o a.out && ./a.out; echo "exit: $?"
