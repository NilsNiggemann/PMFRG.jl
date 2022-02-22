"""
includes MultiLoop corrections to PMFRG.
"""
module MultiLoopPMFRG
    using SpinFRGLattices,OrdinaryDiffEq,DiffEqCallbacks,Parameters,Printf,RecursiveArrayTools,LoopVectorization

    # export Multiloop,SolveFRG,SolveFRG_Checkpoint,generateFileName
    export Multiloop,SetupParquet

    using ..PMFRG: Params,Workspace_Struct,setZero!,iG_,V_,getDFint!,get_Self_Energy!,getVertexDeriv!,mixedFrequencies,double,CreateX,Observables,writeOutput,bufferV_!,getChi,launchPMFRG!,launchPMFRG_Checkpoint,InitializeState,readState,getFileParams,PMFRGMethod
    #functions that are extended
    import ..PMFRG: SolveFRG,SolveFRG_Checkpoint,generateFileName

    struct Multiloop <: PMFRGMethod 
        l::Int
    end

    generateFileName(Par::Params,Method::Multiloop,arg::String = "") = generateFileName(Par,string("_",Method.l,arg))

    struct ParquetWorkspace{T}
        fint::Array{T,1}
        gamma::Array{T,2}
        
        Gamma_a::Array{T,4}
        Gamma_b::Array{T,4}
        Gamma_c::Array{T,4}
        
        Xa::Array{T,4}
        Xb::Array{T,4}
        Xc::Array{T,4}
        
        XTa::Array{T,4}
        XTb::Array{T,4}
        XTc::Array{T,4}
        XTd::Array{T,4}
        
        B0::Array{T,4}
        BX::Array{T,4}
    end
    
    function SetupParquet(Par)
        @unpack N, VDims, System,Ngamma = Par
        @unpack NUnique,Npairs = System
        Workspace = ParquetWorkspace(
            zeros(NUnique), # fint
            zeros(NUnique,Ngamma), # gamma
            zeros(VDims), # Gamma_a
            zeros(VDims), # Gamma_b
            zeros(VDims), # Gamma_c
            
            zeros(VDims), # Xa
            zeros(VDims), # Xb
            zeros(VDims), # Xc
            
            zeros(VDims), # XTa
            zeros(VDims), # XTb
            zeros(VDims), # XTc
            zeros(VDims), # XTd
            
            zeros(VDims), # B0
            zeros(VDims) # BX
        )

    end
    
    """Given set of parameters solve the self-consistent Parquet equations iteratively."""
    function SolveParquet(Par::Params)
    end

    """Obtains a solution to Bethe-Salpeter and Schwinger-Dyson equations by iteration until convergence is reached up to accuracy specified by accuracy in Params"""
    function iterateSolution()

    end

    function iterateSDE(Workspace,Par)

    end
end # module