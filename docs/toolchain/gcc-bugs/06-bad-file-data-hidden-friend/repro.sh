#!/bin/sh
# GCC 16.2: "failed to read compiled module cluster N: Bad file data" /
# "failed to load binding '::f@c'" when importing module c.
# Expected: everything compiles and the program exits with 0.
set -x
CXX=${CXX:-g++}
FLAGS="-std=c++26 -fmodules -Wall -Wextra"
rm -rf gcm.cache
$CXX $FLAGS -c a.cppm -o a.o
$CXX $FLAGS -c b.cppm -o b.o
$CXX $FLAGS -c c.cppm -o c.o
$CXX $FLAGS -c main.cpp -o main.o   # fails here
$CXX a.o b.o c.o main.o -o repro && ./repro
