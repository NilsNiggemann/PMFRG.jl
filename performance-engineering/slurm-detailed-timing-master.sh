#!/bin/bash

submit() {
    local NCPUS=$1
    sbatch --cpus-per-task $NCPUS \
           --output detailed-timing-$NCPUS.out \
           PMFRG.jl/performance-engineering/slurm-benchmarking.sh
}


#for NCPUS in 3 5 6 7 10 15 25 30 50 60 100 120 150 
for NCPUS in 1 2 4 8 19 38 76 152
do
	submit $NCPUS
done
