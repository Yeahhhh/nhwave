#!/bin/env sh

cd ../examples/submerged_bar
time mpirun -np 1 ../../source/build/psg-pgi-debug/nhwave

