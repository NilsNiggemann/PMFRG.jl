module PMFRG
    using SpinFRGLattices,OrdinaryDiffEq,DiffEqCallbacks,Parameters,Printf,RecursiveArrayTools,LoopVectorization,StructArrays
    export SolveFRG,Params,double,BS3,Vern7,DP5,version,getChi

    version() = """v.1.0.7"""
    # Essentials
    include("Globals.jl")
    include("VertexFunctions_Dense.jl")
    include("Propagators.jl")
    include("Flowequations_Dense.jl")
    include("Solver.jl")
    #Precompilation
    # include("precompile.jl")
    
end # module
##
