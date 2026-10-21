#!/bin/sh
# GCC 16.2: undefined reference to __tu_has_violation from an imported inline
# function with a precondition. Expected: links and exits 0.
set -x
CXX=${CXX:-g++}
FLAGS="-std=c++26 -fmodules -Wall -Wextra"
rm -rf gcm.cache *.o a.out

$CXX $FLAGS -c m.cppm -o m.o
$CXX $FLAGS -c main.cpp -o main.o
$CXX m.o main.o -lstdc++exp -o a.out && ./a.out; echo "exit: $?"

# Same failure at -O1 and -O2, and also with
#   -fcontract-disable-optimized-checks -fcontract-checks-outlined
# (which in addition warns "'optimize' attribute currently unsupported in modules"):
for EXTRA in "-O2" "-O2 -fcontract-disable-optimized-checks -fcontract-checks-outlined"; do
  rm -rf gcm.cache *.o a.out
  $CXX $FLAGS $EXTRA -c m.cppm -o m.o
  $CXX $FLAGS $EXTRA -c main.cpp -o main.o
  $CXX m.o main.o -lstdc++exp -o a.out && ./a.out; echo "exit: $?"
done
