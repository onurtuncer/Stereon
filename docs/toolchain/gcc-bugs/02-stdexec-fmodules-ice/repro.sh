#!/bin/sh
# GCC 16.1/16.2: ICE (segfault) in transfer_defining_module; no module code involved.
# Duplicate of PR c++/126783, fixed on trunk and releases/gcc-16 on 2026-08-13 (in 16.3).
set -x
g++ -std=c++26 -fmodules -Wall -Wextra -c main.cpp -o main.o   # ICE on 16.1.0, 16.2.0
g++ -std=c++26           -Wall -Wextra -c main.cpp -o main.o   # OK: needs -fmodules
