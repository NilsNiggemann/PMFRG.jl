module PMFRG
    using SpinFRGLattices,OrdinaryDiffEq,DiffEqCallbacks,Parameters,Printf,RecursiveArrayTools,LoopVectorization,StructArrays,HDF5,H5Zblosc

    using SpinFRGLattices.StaticArrays

    export SolveFRG,Params,OneLoopParams,double,BS3,Vern7,DP5,version,getChi,OneLoop
    version() = """v.2.1.2"""
    # Essentials
    include("Types.jl")
    include("OneLoopTypes.jl")
    include("VertexUtils.jl")

    include("VertexFunctions_Dense.jl")
    include("Propagators.jl")
    include("Flowequations_Dense.jl")
    include("FileIO.jl")
    include("Solver.jl")
    export saveState, readState, readLam, saveParams, readParams, setupDirectory, saveCurrentState, UniqueDirName, generateName, setupFromCheckpoint, SolveFRG_Checkpoint,readGeometry, readObservables,getUnfinishedJobs,generateFileName
    
    
    include("TwoLoop/TwoLoopPMFRG.jl")
    # using .TwoLoopPMFRG
    export TwoLoop
    
    using FixedPoint
    
    include("MultiLoop/MultiLoopPMFRG.jl")
    export MultiLoop,Parquet,SolveParquet

    include("UnitTests/UnitTests.jl")
    # export UnitTests

    #Precompilation
    # include("precompile.jl")
    # __precompile__quiet__()
end # module
##