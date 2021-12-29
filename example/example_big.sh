#!/usr/bin/env bash

mkdir -p results_big
cd results_big
../../traj_analyzer ../traj/CDP_224_big.xyz ../traj/pbc_CDP_224 ../calc.inp
