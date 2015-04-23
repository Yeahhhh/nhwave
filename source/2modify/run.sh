#!/bin/env sh

make clean-all

for config in          \
    psg-gnu-debug      \
    psg-gnu-optimize   \
    psg-intel-debug    \
    psg-intel-optimize \
    psg-pgi-debug      \
    psg-pgi-optimize
do
    make CONFIG="$config" test
done

