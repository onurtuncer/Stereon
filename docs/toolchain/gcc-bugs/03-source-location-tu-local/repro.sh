#!/bin/sh
# Expected (the reported failure): "exposes TU-local entity 'struct std::source_location'"
# while compiling m.cppm. Not reproduced on Compiler Explorer (Linux, GCC 16.2):
# this testcase builds and runs there.
set -e
CXX=${CXX:-g++}
FLAGS="-std=c++26 -fmodules -Wall -Wextra"
$CXX $FLAGS -c m.cppm -o m.o
$CXX $FLAGS -c main.cpp -o main.o
$CXX m.o main.o -o a.out
./a.out
