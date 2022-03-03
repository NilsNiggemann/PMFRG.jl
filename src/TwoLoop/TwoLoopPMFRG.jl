"""
includes Two-loop corrections to PMFRG.
"""
module TwoLoopPMFRG

    using SpinFRGLattices,OrdinaryDiffEq,DiffEqCallbacks,Parameters,Printf,RecursiveArrayTools,LoopVectorization
    export TwoLoop,SolveFRG,SolveFRG_Checkpoint,generateFileName

    using ..PMFRG: Params,Workspace_Struct,setZero!,iG_,V_,getDFint!,get_Self_Energy!,getVertexDeriv!,symmetrizeX!,mixedFrequencies,double,CreateX,CreateXT,Observables,writeOutput,bufferV_!,getChi,launchPMFRG!,launchPMFRG_Checkpoint,InitializeState,readState,getFileParams,PMFRGMethod
    #functions that are extended
    import ..PMFRG: SolveFRG,SolveFRG_Checkpoint,generateFileName,TwoLoop

    function AllocateSetup(Par::Params)
        @unpack N,Ngamma,Npairs,VDims,couplings,T,NUnique = Par
        println("TwoLoop: T= ",T)
        X = CreateX(VDims)
        XTilde = CreateXT(VDims)
        
        Y = CreateX(VDims)
        YTilde = CreateXT(VDims)
        PropsBuffers = [Matrix{double}(undef,NUnique,NUnique) for _ in 1:Threads.nthreads()] 
        VertexBuffers = [VertexBuffer(Par.Npairs) for _ in 1:Threads.nthreads()] 
    
        return (X,XTilde,Y,YTilde,PropsBuffers,VertexBuffers,Par)
    end
    #Overwrite getDeriv function
    include("Flowequations.jl")
    """
    Solves FRG as specified for parameters
    """
    SolveFRG(Par::Params,Method::TwoLoop;kwargs...) = launchPMFRG!(InitializeState(Par),AllocateSetup(Par),getDeriv!; kwargs...)
    SolveFRG_Checkpoint(Filename::String,Geometry::SpinFRGLattices.Geometry,Method::TwoLoop,Par=nothing;kwargs...)= launchPMFRG_Checkpoint(Filename,Geometry,AllocateSetup,getDeriv!,Par;kwargs...)
    generateFileName(Par::Params,Method::TwoLoop,arg::String = "") = generateFileName(Par,"_l2"*arg)
end # module