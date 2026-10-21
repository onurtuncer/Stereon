#!/bin/sh
# ICE in mangle_module when an importer uses a module's defaulted hidden-friend operator==.
set -x
g++ -v
g++ -std=c++20 -fmodules -c m.cppm
g++ -std=c++20 -fmodules -c main.cpp   # ICE here
g++ m.o main.o -o a.out && ./a.out
