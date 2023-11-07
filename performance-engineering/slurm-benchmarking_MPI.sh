#!/bin/bash
#=
#SBATCH --partition cpuonly
#SBATCH --time 5 
#SBATCH --nodes 2
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task=76
#SBATCH --exclusive
#SBATCH --dependency singleton
#SBATCH --job-name pmfrg-benchmark

ROOT="/home/hk-project-scs/hs2454/PMFRG"
PROJECT="$ROOT/TestProject"

set -o nounset
module use "$HOME/modules"
# This is needed anyway by mpiexecjl,
# and needs to match the content of $PROJECT/LocalPreferences.toml
# (as set by MPIPreferences.jl).
module load mpi/openmpi/4.1 
module load julia/1.9.3

MPIEXEC="/home/hk-project-scs/hs2454/.julia/bin/mpiexecjl --project=$PROJECT"
SCRIPT="$ROOT/PMFRG.jl/performance-engineering/slurm-benchmarking_MPI.sh"

COMMAND=($MPIEXEC -n $SLURM_NTASKS 
         julia # --project="$PROJECT" 
         --optimize=3 
         --threads $SLURM_CPUS_PER_TASK 
	 $SCRIPT) 
echo ${COMMAND[@]} 
${COMMAND[@]} 

wait
exit

=#
using MPI
MPI.Init()

rank = 0
nranks = 1

if MPI.Initialized()
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    nranks = MPI.Comm_size(MPI.COMM_WORLD)
end

macro mpi_synchronize(expr)
    quote
    for r in 0:(nranks-1)
        if rank == r
            print("[$rank/$nranks]: ")
            $(esc(expr))
        end
        if MPI.Initialized()
            MPI.Barrier(MPI.COMM_WORLD)
        end
    end
    end
end

function print_barrier(args...)
    @mpi_synchronize println(args...)
end

workdir = "dir$rank-$(Threads.nthreads())"
print_barrier("Removing data from previous runs ($workdir)")
rm(workdir, recursive=true, force=true) 
mkdir(workdir)
cd(workdir)

using ThreadPinning
pinthreads(:cores)


print_barrier("Loading SpinFRGLattices")
using SpinFRGLattices
print_barrier("Loading PMFRG")
using PMFRG
using SpinFRGLattices.SquareLattice
print_barrier("Loading TimerOutputs")
using TimerOutputs

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

print_barrier("GetSquareLattice - system toy")
SystemToy = getSquareLattice(NLenToy, couplings)

print_barrier("GetSquareLattice")
System = getSquareLattice(NLen, couplings) 

print_barrier("Warm up")

print_barrier("Get Params - toy")
Par = Params( #create a group of all parameters to pass them to the FRG Solver
    SystemToy, # geometry, this is always required
    OneLoop(), # method. OneLoop() is the default
    T=0.5, # Temperature for the simulation.
    N=10, # Number of positive Matsubara frequencies for the four-point vertex.
    accuracy=1e-3, #absolute and relative tolerance of the ODE solver.
    # For further optional arguments, see documentation of 'NumericalParams'
    MinimalOutput=true,
)

tempdir = "temp-$rank"
print_barrier("Removing data from previous runs ($tempdir)")
rm(tempdir, recursive=true, force=true)
mainFile = "$tempdir/" * PMFRG.generateFileName(Par, "_testFile") # specify a file name for main Output
flowpath = "$tempdir/flows/" # specify path for vertex checkpoints

print_barrier("SolveFRG - toy")
Solution, saved_values = SolveFRG(
    Par,
    MainFile=mainFile,
    CheckpointDirectory=flowpath,
    method=DP5(),
    VertexCheckpoints=[],
    CheckPointSteps=3,
);



print_barrier("Warmup done, timing real problem now.")


print_barrier("Get Params")
Par = Params( #create a group of all parameters to pass them to the FRG Solver
    System, # geometry, this is always required
    OneLoop(), # method. OneLoop() is the default
    T=0.5, # Temperature for the simulation.
    N=25, # Number of positive Matsubara frequencies for the four-point vertex.
    accuracy=1e-3, #absolute and relative tolerance of the ODE solver.
    # For further optional arguments, see documentation of 'NumericalParams'
    MinimalOutput=true,
)

tempdir = "temp-$rank"
print_barrier("Removing data from previous runs ($tempdir)")
rm(tempdir, recursive=true, force=true)
mainFile = "$tempdir/" * PMFRG.generateFileName(Par, "_testFile") # specify a file name for main Output
flowpath = "$tempdir/flows/" # specify path for vertex checkpoints

reset_timer!()
print_barrier("SolveFRG")
@time Solution, saved_values = SolveFRG(
    Par,
    MainFile=mainFile,
    CheckpointDirectory=flowpath,
    method=DP5(),
    VertexCheckpoints=[],
    CheckPointSteps=3,
);

@mpi_synchronize print_timer()

if MPI.Initialized()
 MPI.Finalize()
end



