"""
includes Two-loop corrections to PMFRG.
"""
# module TwoLoopPMFRG

#     using SpinFRGLattices,OrdinaryDiffEq,DiffEqCallbacks,Parameters,Printf,RecursiveArrayTools,LoopVectorization
#     export TwoLoop,SolveFRG,SolveFRG_Checkpoint,generateFileName,AllocateSetup

#     using ..PMFRG: BubbleType,StateType,VertexType,PMFRGWorkspace,PMFRGParams,OneLoopWorkspace,NumericalParams,OptionParams,Observables,double
    
#     using ..PMFRG: setZero!,iG_,V_,getDFint!,get_Self_Energy!,getVertexDeriv!,symmetrizeBubble!,mixedFrequencies,writeOutput,bufferV_!,getChi,launchPMFRG!,launchPMFRG_Checkpoint,InitializeState,readState,getFileParams
    
#     #functions that are extended
#     import ..PMFRG: SolveFRG,SolveFRG_Checkpoint,generateFileName,Params,AllocateSetup

    include("TwoLoopTypes.jl")

    function AllocateSetup(Par::TwoLoopParams)
        @unpack N,Ngamma = Par.NumericalParams
    
        @unpack couplings,NUnique,Npairs = Par.System
        println("Two Loop: T= ",Par.NumericalParams.T)
        ##Allocate Memory:
        X = BubbleType(Par)
        Y = BubbleType(Par)
        PropsBuffers = [Matrix{double}(undef,NUnique,NUnique) for _ in 1:Threads.nthreads()] 
        VertexBuffers = [VertexBufferType(Npairs) for _ in 1:Threads.nthreads()]
        BubbleBuffers = [BubbleBufferType(Npairs) for _ in 1:Threads.nthreads()]
        Buffs = BufferTypeTwoLoop(PropsBuffers,VertexBuffers,BubbleBuffers) 
        return (X,Y,Buffs,Par)
    end
    include("Flowequations.jl")

    Params(System::Geometry,O::TwoLoop;kwargs...) = TwoLoopParams(System,NumericalParams(;kwargs...),OptionParams(;kwargs...))

    """
    Solves FRG as specified for parameters
    """
    SolveFRG(Par::TwoLoopParams;kwargs...) = launchPMFRG!(InitializeState(Par),AllocateSetup(Par),getDeriv2L!; kwargs...)

    SolveFRG_Checkpoint(Filename::String,Geometry::SpinFRGLattices.Geometry,Method::TwoLoop,Par=nothing;kwargs...)= launchPMFRG_Checkpoint(Filename,Geometry,AllocateSetup,getDeriv2L!,Par;kwargs...)

    generateFileName(Par::TwoLoopParams,Method::TwoLoop,arg::String = "") = generateFileName(Par,"_l2"*arg)
# end # module