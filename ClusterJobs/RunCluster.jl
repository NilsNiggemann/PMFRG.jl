include("ClusterUtils.jl")
# include("../src/PMFRG.jl")
@time @everywhere using PMFRG
using SpinFRGLattices
# include("../TwoLoopPMFRG.jl")
# using .SimpleCubic
##
System = getPolymer(4)
DistributedTSweep(0.5:0.1:2,System,N=24,Ngamma = 300,add_name = "2",accuracy = 1E-9)