#!/bin/bash
#=
#SBATCH --partition cpuonly
#SBATCH --time 25
#SBATCH --nodes 1
#SBATCH --ntasks 2
#SBATCH --cpus-per-task=38
#SBATCH --exclusive
#SBATCH --dependency singleton
#SBATCH --job-name pmfrg-benchmark

set -o nounset
module use "$HOME/modules"
module load julia/1.9.3
module load mpi/openmpi/4.1

rm -rf  dir$SLURM_CPUS_PER_TASK 

# From this discussion:
# https://discourse.julialang.org/t/compilation-options-for-downfall-mitigation/104844
# the -Cnative,-fast-gather options should help 
# with the downfall mitigation issue for the moment
mpirun -n $SLURM_NTASKS \
      julia --optimize=3 \
      --cpu-target native,-fast-gather \
      --threads $SLURM_CPUS_PER_TASK \
      ${BASH_SOURCE[0]} &

wait
exit

=#
import Pkg
ROOT = "/home/hk-project-scs/hs2454/PMFRG/"
Pkg.activate(ROOT * "TestProject" )
using SpinFRGLattices, PMFRG
using SpinFRGLattices.SquareLattice
using TimerOutputs
using MPI
using ThreadPinning

MPI.Init()
pinthreads(:cores)
threadinfo(color=false)

rank = MPI.Comm_rank(MPI_COMM_WORLD)

workdir = "dir$rank-$(Threads.nthreads())"
println("Removing data from previous runs ($workdir)")
rm(orkdir, recursive=true, force=true) 
mkdir(workdir)
cd(workdir)
TimerOutputs.enable_debug_timings(PMFRG)

# Number of nearest neighbor bonds 
# up to which correlations are treated in the lattice. 
# For NLen = 5, all correlations C_{ij} are zero 
#if sites i and j are separated by more than 5 nearest neighbor bonds.
NLenToy = 5 
NLen = 14 
J1 = 1
J2 = 0.1
# Construct a vector of couplings: 
# nearest neighbor coupling is J1 (J2) 
# and further couplings to zero.
# For finite further couplings simply provide a longer array, 
# i.e [J1,J2,J3,...]
couplings = [J1, J2] 

# create a structure that contains all information about the geometry of the problem.

SystemToy = getSquareLattice(NLenToy, couplings)

System = getSquareLattice(NLen, couplings) 

println("Warm up")

Par = Params( #create a group of all parameters to pass them to the FRG Solver
    SystemToy, # geometry, this is always required
    OneLoop(), # method. OneLoop() is the default
    T=0.5, # Temperature for the simulation.
    N=10, # Number of positive Matsubara frequencies for the four-point vertex.
    accuracy=1e-3, #absolute and relative tolerance of the ODE solver.
    # For further optional arguments, see documentation of 'NumericalParams'
    MinimalOutput=true,
)

tempdir = "temp"
mainFile = "$tempdir/" * PMFRG.generateFileName(Par, "_testFile") # specify a file name for main Output
flowpath = "$tempdir/flows/" # specify path for vertex checkpoints

Solution, saved_values = SolveFRG(
    Par,
    MainFile=mainFile,
    CheckpointDirectory=flowpath,
    method=DP5(),
    VertexCheckpoints=[],
    CheckPointSteps=3,
);



println("Warmup done, timing real problem now.")


Par = Params( #create a group of all parameters to pass them to the FRG Solver
    System, # geometry, this is always required
    OneLoop(), # method. OneLoop() is the default
    T=0.5, # Temperature for the simulation.
    N=25, # Number of positive Matsubara frequencies for the four-point vertex.
    accuracy=1e-3, #absolute and relative tolerance of the ODE solver.
    # For further optional arguments, see documentation of 'NumericalParams'
    MinimalOutput=true,
)

tempdir = "temp"
println("Removing data from previous runs ($tempdir)")
rm(tempdir, recursive=true, force=true)
mainFile = "$tempdir/" * PMFRG.generateFileName(Par, "_testFile") # specify a file name for main Output
flowpath = "$tempdir/flows/" # specify path for vertex checkpoints

reset_timer!()
@time Solution, saved_values = SolveFRG(
    Par,
    MainFile=mainFile,
    CheckpointDirectory=flowpath,
    method=DP5(),
    VertexCheckpoints=[],
    CheckPointSteps=3,
);

print_timer()
println("Done")
MPI.Finalize()

println("Pinning after the run:")
threadinfo(color=false)


