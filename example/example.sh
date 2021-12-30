#!/usr/bin/env bash

mkdir -p results
cd results
../../traj_analyzer ../traj/CDP_224.xyz ../traj/pbc_CDP_224 ../calc.inp
#valgrind --leak-check=full ../../traj_analyzer ../traj/CDP_224.xyz ../traj/pbc_CDP_224 ../calc.inp
