"""
includes MultiLoop corrections to PMFRG.
"""
module MultiLoopPMFRG
    using SpinFRGLattices,OrdinaryDiffEq,DiffEqCallbacks,Parameters,Printf,RecursiveArrayTools,LoopVectorization,StructArrays

    # export MultiLoop,SolveParquet

    using ..PMFRG: MultiLoop,Params,Workspace_Struct,iG_,V_,getDFint!,get_Self_Energy!,getXBubble!,mixedFrequencies,double,CreateX,Observables,bufferV_!,convertFreqArgs,getChi,launchPMFRG!,launchPMFRG_Checkpoint,InitializeState,readState,getFileParams,PMFRGMethod
    #functions that are extended
    import ..PMFRG: SolveFRG,SolveFRG_Checkpoint,generateFileName,getChi,writeOutput
    include("Types.jl")
    include("VertexUtils.jl")
    include("Flowequations.jl")
    include("Parquet.jl")

    export SolveParquet
    generateFileName(Par::Params,Method::MultiLoop,arg::String = "") = generateFileName(Par,string("_",Method.l,arg))

end # module