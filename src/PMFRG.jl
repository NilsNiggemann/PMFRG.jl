module PMFRG
    using SpinFRGLattices,OrdinaryDiffEq,DiffEqCallbacks,Parameters,Printf,RecursiveArrayTools,LoopVectorization
    export SolveFRG,Params,double,BS3,Vern7,DP5,version

    version() = """v.1.0.2"""
    # Essentials
    include("Globals.jl")
    include("VertexFunctions_Dense.jl")
    include("Propagators.jl")
    include("Flowequations_Dense.jl")
    include("Solver.jl")
    
    function Evaluation_Mode() 
        @eval Main begin

            include(string($(@__DIR__),"/Evaluation.jl"))
        end
    end
    #Precompilation
    # include("precompile.jl")
    
end # module
##
