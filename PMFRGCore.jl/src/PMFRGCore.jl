module PMFRGCore
using SpinFRGLattices, Printf, LoopVectorization, StructArrays, HDF5, H5Zblosc, TimerOutputs

using SpinFRGLattices.StaticArrays

export Params, OneLoopParams, version, getChi, OneLoop
version() = v"2.2.1"
# Essentials
include("Types.jl")
include("OneLoopTypes.jl")
include("VertexUtils.jl")
include("StateLib.jl")
using .StateLib
include("VertexFunctions_Dense.jl")
include("Propagators.jl")
include("Flowequations_Dense.jl")
include("FileIO.jl")
include("SolverUtils.jl")
export saveState,
    readState,
    readLam,
    saveParams,
    readParams,
    setupDirectory,
    UniqueDirName,
    UniqueFileName,
    generateName,
    readGeometry,
    readObservables,
    getUnfinishedJobs,
    generateFileName,
    generateMainFile,
    AbstractParallelizationScheme,
    MultiThreaded,
    PMFRGParams

include("TwoLoop/TwoLoopPMFRG.jl")
# using .TwoLoopPMFRG
export TwoLoop

using FixedPoint

include("MultiLoop/MultiLoopPMFRG.jl")
export MultiLoop, Parquet

export UseMPI

# export UnitTests

#Precompilation
# include("precompile.jl")
# __precompile__quiet__()
end # module
##
