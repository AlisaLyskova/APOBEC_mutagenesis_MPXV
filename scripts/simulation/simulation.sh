#!/bin/bash

#for DNA samples
OUT_SHARES=simulations_shares.csv
OUT_NREADS=simulations_Nreads_positions.csv
OUT_LOG=simulations.log
touch $OUT_SHARES
touch $OUT_NREADS
touch $OUT_LOG

sample_arr=( $(python3 features_samples.py) )

parallel -j 15 python3 simulation.py {1} ::: "${sample_arr[@]}"


#for dRNA samples
sample_arr=( $(python3 features_samples_dRNA.py) )

parallel -j 10 python3 simulation.py {1} ::: "${sample_arr[@]}"
