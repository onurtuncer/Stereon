#!/bin/sh
# MSYS2 UCRT64 GCC 16.2 (x86_64-w64-mingw32): link fails with
# "multiple definition of <function-local static of a libstdc++ inline function>".
# Linux: links, but the copy emitted in a.o is WEAK instead of UNIQUE (see readelf).
set -x
rm -rf gcm.cache *.o *.s a.out a.exe
g++ -std=c++26 -fmodules -Wall -Wextra -x c++ -c m.cppm -o m.o
g++ -std=c++26 -fmodules -Wall -Wextra -c a.cpp -o a.o
g++ -std=c++26 -fmodules -Wall -Wextra -c main.cpp -o main.o
g++ m.o a.o main.o          # MinGW: multiple definition ... first defined here
# Root-cause evidence: how each TU emits the static.
g++ -std=c++26 -fmodules -S a.cpp -o a.s
g++ -std=c++26 -fmodules -S main.cpp -o main.s
grep -n -B2 -A3 -E '(__replace_rep|_S_tiEvE5__tag)"?:' a.s main.s
#   MinGW: expect main.s ".section .rdata$...,\"dr\"" + ".linkonce same_size"; a.s without .linkonce?
#   ELF:   expect main.s "@gnu_unique_object" + .comdat; a.s ".weak" only?
readelf -sW a.o main.o 2>/dev/null | grep -E '__replace_rep|_S_tiEvE5__tag'
