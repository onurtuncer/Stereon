#!/bin/sh
# Runs every repro.sh in a scratch copy and keeps each log next to it as repro.log.
# Usage: sh run-all.sh   (from docs/toolchain/gcc-bugs, with GCC 16 on PATH)
here=$(cd "$(dirname "$0")" && pwd)
g++ --version | head -n 1
for d in "$here"/[0-9][0-9]-*/; do
    [ -f "$d/repro.sh" ] || continue
    name=$(basename "$d")
    work=$(mktemp -d)
    cp "$d"/* "$work"/
    (cd "$work" && sh repro.sh) > "$d/repro.log" 2>&1
    echo "$name: exit $? (log: $name/repro.log)"
    rm -rf "$work"
done
