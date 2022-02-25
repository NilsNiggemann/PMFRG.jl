"""
includes MultiLoop corrections to PMFRG.
"""
module MultiLoopPMFRG
    using SpinFRGLattices,OrdinaryDiffEq,DiffEqCallbacks,Parameters,Printf,RecursiveArrayTools,LoopVectorization

    # export Multiloop,SolveFRG,SolveFRG_Checkpoint,generateFileName
    export Multiloop,SetupParquet

    using ..PMFRG: Multiloop,State,Bubble,Params,Workspace_Struct,setZero!,iG_,V_,getDFint!,get_Self_Energy!,getVertexDeriv!,mixedFrequencies,double,CreateX,Observables,writeOutput,bufferV_!,getChi,launchPMFRG!,launchPMFRG_Checkpoint,InitializeState,readState,getFileParams,PMFRGMethod
    #functions that are extended
    import ..PMFRG: SolveFRG,SolveFRG_Checkpoint,generateFileName

    generateFileName(Par::Params,Method::Multiloop,arg::String = "") = generateFileName(Par,string("_",Method.l,arg))

    struct ParquetWorkspace{T,P <: Params{S,M} where {S,M}}
        State::StateType{T} #Stores the current step of iteration
        NewState::StateType{T} #Stores the next step of iteration
        
        X::BubbleType{T} #Stores the actual bubble funtion X = B0 + BX

        B0::BubbleType{T} # Bubble involving bare vertex Γ_0 ∘ Γ
        BX::BubbleType{T} # Bubble involving other bubble X ∘ Γ
        Par::P # Params
    end
    
    function SetupParquet(Par)
        Workspace = ParquetWorkspace(
            State(Par),
            State(Par),
            Bubble(Par),

            Bubble(Par),
            Bubble(Par),

            Par
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