#!/bin/sh
# GCC 16.2: pre/post of a function template imported from a module are not
# checked. Expected: contract violation reported, program aborts (exit 134).
# Actual: exits 0.
set -x
CXX=${CXX:-g++}
FLAGS="-std=c++26 -fmodules -Wall -Wextra"
rm -rf gcm.cache *.o a.out

$CXX $FLAGS -c m.cppm -o m.o
$CXX $FLAGS -c main.cpp -o main.o
$CXX m.o main.o -lstdc++exp -o a.out
./a.out; echo "exit: $? (expected 134, a contract violation)"

# Control: the same template without modules is checked.
cat > control.cpp <<'EOF'
template <class T> T g(T x) pre(x > 0) { return x; }
int main(int argc, char **) { return g(argc - 1); }
EOF
$CXX -std=c++26 -Wall -Wextra control.cpp -lstdc++exp -o control
./control; echo "control exit: $?"
