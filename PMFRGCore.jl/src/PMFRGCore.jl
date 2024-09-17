module PMFRGCore
using SpinFRGLattices,
    Printf,
    RecursiveArrayTools,
    LoopVectorization,
    StructArrays,
    HDF5,
    H5Zblosc,
    SciMLBase,
    DiffEqCallbacks,
    TimerOutputs

using SpinFRGLattices.StaticArrays

export Params, OneLoopParams, version, getChi, OneLoop
version() = v"2.2.1"
# Essentials
include("Types.jl")
include("OneLoopTypes.jl")
include("VertexUtils.jl")

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
    saveCurrentState,
    UniqueDirName,
    UniqueFileName,
    generateName,
    readGeometry,
    readObservables,
    getUnfinishedJobs,
    generateFileName,
    generateMainFile,
    AbstractParallelizationScheme,
    MultiThreaded

include("TwoLoop/TwoLoopPMFRG.jl")
# using .TwoLoopPMFRG
export TwoLoop

using FixedPoint

include("MultiLoop/MultiLoopPMFRG.jl")
export MultiLoop, Parquet, SolveParquet

export UseMPI

# export UnitTests

#Precompilation
# include("precompile.jl")
# __precompile__quiet__()
end # module
##
