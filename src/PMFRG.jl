module PMFRG
    using SpinFRGLattices,OrdinaryDiffEq,DiffEqCallbacks,Parameters,Printf,RecursiveArrayTools,LoopVectorization,StructArrays,HDF5
    export SolveFRG,Params,double,BS3,Vern7,DP5,version,getChi,OneLoop
    version() = """v.2.1.0"""
    # Essentials
    include("Globals.jl")
    include("VertexFunctions_Dense.jl")
    include("Propagators.jl")
    include("Flowequations_Dense.jl")
    include("FileIO.jl")
    include("Solver.jl")
    export saveState, readState, readLam, saveParams, readParams, setupDirectory, saveCurrentState, UniqueDirName, generateName, setupFromCheckpoint, SolveFRG_Checkpoint,readGeometry, readObservables,getUnfinishedJobs,generateFileName
    
    include("TwoLoop/TwoLoopPMFRG.jl")
    using .TwoLoopPMFRG
    export TwoLoop
    
    include("MultiLoop/MultiLoopPMFRG.jl")
    using .MultiLoopPMFRG
    export MultiLoop,SolveParquet
    # Precompilation
    # include("precompile.jl")
    # __precompile__quiet__()
end # module
##