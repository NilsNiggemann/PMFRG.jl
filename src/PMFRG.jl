module PMFRG
    using SpinFRGLattices,OrdinaryDiffEq,DiffEqCallbacks,Parameters,Printf,RecursiveArrayTools,OrdinaryDiffEq,DiffEqCallbacks,Parameters,Printf
    export SolveFRG,Params

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
