#!/bin/sh
# GCC: undefined reference to an inline explicit specialization of a
# variable template odr-used by an importer.
# Expected: links and exits 0. Actual: undefined reference to `N@a<L@a>'.
set -x
CXX=${CXX:-g++}
FLAGS="-std=c++26 -fmodules -Wall -Wextra"
rm -rf gcm.cache
$CXX $FLAGS -c a.cppm
$CXX $FLAGS -c main.cpp
$CXX a.o main.o -o main && ./main; echo "exit: $?"
