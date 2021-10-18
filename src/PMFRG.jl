module PMFRG
    using SpinFRGLattices,OrdinaryDiffEq,DiffEqCallbacks,Parameters,Printf,RecursiveArrayTools,LoopVectorization,StructArrays,HDF5
    export SolveFRG,Params,double,BS3,Vern7,DP5,version,getChi

    version() = """v.1.1.0"""
    # Essentials
    include("Globals.jl")
    include("VertexFunctions_Dense.jl")
    include("Propagators.jl")
    include("Flowequations_Dense.jl")
    include("FileIO.jl")
    include("Solver.jl")
    export saveState, readState, readLam, saveParams, readParams, setupDirectory, saveCurrentState, UniqueDirName, generateName, setupFromCheckpoint, SolveFRG_Checkpoint
    #Precompilation
    # include("precompile.jl")
    
end # module
##