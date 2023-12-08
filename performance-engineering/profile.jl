#!/usr/bin/env julia
# This example comes from README.md

using Pkg
using SpinFRGLattices, PMFRG
using SpinFRGLattices.SquareLattice

using Profile

print("All modules loaded")
# Number of nearest neighbor bonds
# up to which correlations are treated in the lattice.
# For NLen = 5, all correlations C_{ij} are zero
# if sites i and j are separated by more
# than 5 nearest neighbor bonds.

NLenToy = 5
NLen = 14

# Construct a vector of couplings: nearest neighbor coupling
# is J1 (J2) and further couplings to zero.
# For finite further couplings simply provide
# a longer array, i.e [J1,J2,J3,...]
J1 = 1
J2 = 0.5
couplings = [J1, J2]

# create a structure that contains all information about the geometry of the problem.
SystemToy = getSquareLattice(NLenToy, couplings)
System = getSquareLattice(NLen, couplings)

#create a group of all parameters to pass them to the FRG Solver
# For further optional arguments, see documentation of 'NumericalParams'
ParSmall = Params(
    SystemToy,        # geometry, this is always required
    OneLoop(),     # method. OneLoop() is the default
    T = 0.5,         # Temperature for the simulation.
    N = 10,          # Number of positive Matsubara frequencies for the four-point vertex.
    accuracy = 1e-3, #absolute and relative tolerance of the ODE solver.
    MinimalOutput = true,
)
Par = Params(
    System,        # geometry, this is always required
    OneLoop(),     # method. OneLoop() is the default
    T = 0.5,         # Temperature for the simulation.
    N = 25,          # Number of positive Matsubara frequencies for the four-point vertex.
    accuracy = 1e-3, #absolute and relative tolerance of the ODE solver.
    MinimalOutput = true,
)



function profile_solvefrg()

    # specify a file name for main Output
    mainFile = "profile-playground/" * PMFRG.generateFileName(ParSmall, "_testFile")

    # specify path for vertex checkpoints
    flowpath = "profile-playground/flows/"

    # First run to ensure we have the compiled code
    # and we are not profiling the jit compiler instead.
    Solution, saved_values = SolveFRG(
        ParSmall,
        MainFile = mainFile,
        CheckpointDirectory = flowpath,
        method = DP5(),
        VertexCheckpoints = [],
        CheckPointSteps = 3,
    )

    rm("profile-playground/", recursive = true)

    mainFile = "profile-playground/" * PMFRG.generateFileName(Par, "_testFile")
    @profile Solution, saved_values = SolveFRG(
        Par,
        MainFile = mainFile,
        CheckpointDirectory = flowpath,
        method = DP5(),
        VertexCheckpoints = [],
        CheckPointSteps = 3,
    )

end


function profile_and_save()
    profile_solvefrg()
    open("profile-threads=$(Threads.nthreads())", "w") do s
        Profile.print(IOContext(s, :displaysize => (24, 500)))
    end
end
