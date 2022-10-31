using SpinFRGLattices, PMFRG
using SpinFRGLattices.SquareLattice

NLen = 5 # Number of nearest neighbor bonds up to which correlations are treated in the lattice. For NLen = 5, all correlations C_{ij} are zero if sites i and j are separated by more than 5 nearest neighbor bonds.
J1 = 1
J2 = 0.1
couplings = [J1,J2] # Construct a vector of couplings: nearest neighbor coupling is J1 (J2) and further couplings to zero. For finite further couplings simply provide a longer array, i.e [J1,J2,J3,...]

System = getSquareLattice(NLen,couplings) # create a structure that contains all information about the geometry of the problem. 

Par = Params( #create a group of all parameters to pass them to the FRG Solver
    System, # geometry, this is always required
    OneLoop(), # method. OneLoop() is the default
    T=0.5, # Temperature for the simulation.
    N = 10, # Number of positive Matsubara frequencies for the four-point vertex.
    accuracy = 1e-3, #absolute and relative tolerance of the ODE solver.
    # For further optional arguments, see documentation of 'NumericalParams'
)

mainFile = "temp/"*PMFRG.generateFileName(Par,"_testFile") # specify a file name for main Output
flowpath = "temp/flows/" # specify path for vertex checkpoints

Solution,saved_values = SolveFRG(Par,MainFile = mainFile ,CheckpointDirectory = flowpath,method = DP5(),VertexCheckpoints = [],CheckPointSteps = 3 )

## Evaluation part

using HDF5, FRGLatticeEvaluation
using CairoMakie #for plotting. You can use whatever plotting package you like of course
##
Lattice = LatticeInfo(System,SquareLattice)

chi_ΛR = h5read(mainFile,"0.5/Chi")
chi_R = chi_ΛR[:,end]
T = h5read(mainFile,"0.5/T")

chi = getFourier(chi_R,Lattice)

k = LinRange(-2pi,2pi,100)

chik = [chi(x,y) for x in k, y in k]

heatmap(k,k,chik)