#!/bin/bash
#=
#SBATCH --partition cpuonly
#SBATCH --time 25
#SBATCH --nodes 1

set -o nounset
module use "$HOME/modules"
module load julia/1.9.2

julia --optimize=3 \
      --threads $SLURM_CPUS_PER_TASK \
      ${BASH_SOURCE[0]}
exit

=#

import Pkg
ROOT = "/home/hk-project-scs/hs2454/PMFRG/"
Pkg.activate(ROOT * "TestProject" )
usin ThreadPinning
println("Before Thread pinning:")
threadinfo()
pinthreads(:core)
println("After Thread pinning:")
threadinfo()

workdir = "dir$(Threads.nthreads())"
mkdir(workdir)
cd(workdir)
include(ROOT * "PMFRG.jl/performance-engineering/generating_data_benchmark.jl")
