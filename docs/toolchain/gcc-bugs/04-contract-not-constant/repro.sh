#!/bin/sh
# GCC 16.1/16.2: "contract condition is not constant" for a precondition that
# calls a member function on an object whose member was already called in the
# same constant evaluation. No modules involved. Fixed on trunk (20261020).
# Expected: compiles, exits 0.
set -x
CXX=${CXX:-g++}
$CXX -std=c++26 -Wall -Wextra main.cpp -lstdc++exp -o a.out && ./a.out; echo "exit: $?"
