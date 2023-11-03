#!/bin/bash
#=
#SBATCH --partition cpuonly
#SBATCH --time 25
#SBATCH --nodes 1
#SBATCH --exclusive
#SBATCH --dependency singleton
#SBATCH --job-name pmfrg-benchmark

set -o nounset
module use "$HOME/modules"
module load julia/1.9.3

rm -rf  dir$SLURM_CPUS_PER_TASK 

# From this discussion:
# https://discourse.julialang.org/t/compilation-options-for-downfall-mitigation/104844
# the -Cnative,-fast-gather options should help 
# with the downfall mitigation issue for the moment
srun julia --optimize=3 \
      --cpu-target native,-fast-gather \
      --threads $SLURM_CPUS_PER_TASK \
      ${BASH_SOURCE[0]} &

wait
exit

=#

import Pkg
ROOT = "/home/hk-project-scs/hs2454/PMFRG/"
Pkg.activate(ROOT * "TestProject" )

using ThreadPinning
pinthreads(:cores)
threadinfo(color=false)

workdir = "dir$(Threads.nthreads())"
mkdir(workdir)
cd(workdir)
include(ROOT * "PMFRG.jl/performance-engineering/generating_data_benchmark.jl")

println("Pinning after the run:")
threadinfo(color=false)


