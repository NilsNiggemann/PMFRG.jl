using SpinFRGLattices, PMFRG
using SpinFRGLattices.SquareLattice

# Number of nearest neighbor bonds
# up to which correlations are treated in the lattice.
# For NLen = 5, all correlations C_{ij} are zero
# if sites i and j are separated by more than 5 nearest neighbor bonds.
NLen = 5
J1 = 1
J2 = 0.1
# Construct a vector of couplings:
# nearest neighbor coupling is J1 (J2)
# and further couplings to zero.
# For finite further couplings
# simply provide a longer array, i.e [J1,J2,J3,...]
couplings = [J1, J2]

# create a structure that contains all information
# about the geometry of the problem.
System = getSquareLattice(NLen, couplings)

#create a group of all parameters to pass them to the FRG Solver
# For further optional arguments, see documentation of 'NumericalParams'
Par = Params(
    System,          # geometry, this is always required
    OneLoop(),       # method. OneLoop() is the default
    T = 0.5,         # Temperature for the simulation.
    N = 10,          # Number of positive Matsubara frequencies for the four-point vertex.
    accuracy = 1e-3, # absolute and relative tolerance of the ODE solver.
)

tempdir = "temp"
println("Removing data from previous runs ($tempdir)")
rm(tempdir, recursive = true, force = true)
# specify a file name for main Output
mainFile = "$tempdir/" * PMFRG.generateFileName(Par, "_testFile")
# specify path for vertex checkpoints
flowpath = "$tempdir/flows/"

Solution, saved_values = SolveFRG(
    Par,
    MainFile = mainFile,
    CheckpointDirectory = flowpath,
    method = DP5(),
    VertexCheckpoints = [],
    CheckPointSteps = 3,
)
