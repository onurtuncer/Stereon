#!/bin/sh
# Wrong code: an importer passes an empty class argument in a stack slot.
# Expected output "123", exit status 0. Fails on x86-64 Linux (SysV ABI).
set -e
CXX=${CXX:-g++}
FLAGS="-std=c++26 -fmodules -Wall -Wextra ${OPT:-}"
$CXX $FLAGS -c m.cppm -o m.o
$CXX $FLAGS -c main.cpp -o main.o
$CXX m.o main.o -o a.out
./a.out
