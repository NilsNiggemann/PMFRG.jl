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

srun julia --optimize=3 \
      --threads $SLURM_CPUS_PER_TASK \
      ${BASH_SOURCE[0]} &

wait
exit

=#

import Pkg
ROOT = "/home/hk-project-scs/hs2454/PMFRG/"
Pkg.activate(ROOT * "TestProject" )
# Printing all dependencies in job log.
deps = Pkg.dependencies()
for dep in deps
    println(dep.first)
    for field in fieldnames(typeof(dep.second))
	values = getfield(dep.second,field)
	print("    ",field)
        if typeof(values) == Dict{String,Base.UUID}
	    for k in values
	        println(k)
            end
        else
            print(values,"\n")		
	end
    end
end

using ThreadPinning
println("Before Thread pinning:")
threadinfo(color=false)
pinthreads(:cores)
println("After Thread pinning:")
threadinfo(color=false)

workdir = "dir$(Threads.nthreads())"
mkdir(workdir)
cd(workdir)
include(ROOT * "PMFRG.jl/performance-engineering/generating_data_benchmark.jl")

println("Pinning after the run:")
threadinfo(color=false)


