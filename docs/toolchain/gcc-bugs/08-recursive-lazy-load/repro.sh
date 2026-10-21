#!/bin/sh
# GCC "recursive lazy load" / "failed to load pendings for '::V'".
# Expected: everything compiles. Actual: c.cppm fails.
set -x
CXX=${CXX:-g++}
FLAGS="-std=c++26 -fmodules -Wall -Wextra"
rm -rf gcm.cache
$CXX $FLAGS -c a.cppm
$CXX $FLAGS -c b.cppm
$CXX $FLAGS -c c-part.cppm
$CXX $FLAGS -c c.cppm      # error: recursive lazy load
$CXX $FLAGS -c main.cpp
$CXX a.o b.o c-part.o c.o main.o -o main && ./main
