module PMFRG
    using SpinFRGLattices,OrdinaryDiffEq,DiffEqCallbacks,Parameters,Printf,RecursiveArrayTools
    export SolveFRG,Params,double,BS3,Vern7,DP5,version

    version() = """v.0.1.1"""
    # Essentials
    include("Globals.jl")
    include("VertexFunctions_Dense.jl")
    include("Propagators.jl")
    include("Flowequations_Dense.jl")
    include("Solver.jl")
    #Precompilation
    include("precompile.jl")
    
    precompile(Testrun(),())
end # module
##
