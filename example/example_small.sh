#!/usr/bin/env bash

mkdir -p results_small
cd results_small
../../dostuff.out ../traj/CDP_224_small.xyz ../traj/pbc_CDP_224 ../calc.inp
