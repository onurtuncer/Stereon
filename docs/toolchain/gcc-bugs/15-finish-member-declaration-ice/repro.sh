#!/bin/sh
# GCC ICE in finish_member_declaration (cp/semantics.cc:4261 in 16.2.0).
# Expected: builds and exits 0. Actual: ICE compiling main.cpp.
set -x
CXX=${CXX:-g++}
FLAGS="-std=c++26 -fmodules -Wall -Wextra"
rm -rf gcm.cache
$CXX $FLAGS -c a.cppm
$CXX $FLAGS -c b.cppm
$CXX $FLAGS -c main.cpp    # internal compiler error
$CXX a.o b.o main.o -o main && ./main; echo "exit: $?"
