struct MultiLoop
    l::Int
end

struct Parquet
end

struct MultiLoopParams{F,G <: Geometry} <: PMFRGParams
    System::G
    NumericalParams::NumericalParams{F}
    Options::OptionParams
    l::Int
end

struct ParquetOptions{F} <: AbstractOptions
    maxIterBSE::Int
    maxIterSDE::Int
    SDE_tolerance::F
    usesymmetry::Bool
    MinimalOutput::Bool
end

function ParquetOptions(;
    maxIterBSE::Int = 40,
    maxIterSDE::Int = 1000,
    SDE_tolerance::AbstractFloat = 1e-9,
    usesymmetry::Bool = true,
    MinimalOutput::Bool = false,
    kwargs...) 
    return ParquetOptions(maxIterBSE,maxIterSDE,SDE_tolerance,usesymmetry,MinimalOutput)
end

struct ParquetParams{F,G <: Geometry} <: PMFRGParams
    System::G
    NumericalParams::NumericalParams{F}
    Options::ParquetOptions{F}
end

"""Todo: make this error when an unknown kwarg is given!"""
Params(System::Geometry,O::MultiLoop;kwargs...) = MultiLoopParams(System,NumericalParams(;kwargs...),OptionParams(;kwargs...),O.l)

Params(System::Geometry,O::Parquet;kwargs...) = ParquetParams(System,NumericalParams(;kwargs...),ParquetOptions(;kwargs...))


getLoopOrder(P::MultiLoopParams) = P.l
getPMFRGMethod(P::MultiLoopParams) = MultiLoop(getLoopOrder(P))
getPMFRGMethod(P::ParquetParams) = Parquet()
generateFileName(Par::MultiLoopParams,arg::String = "") = _generateFileName(Par,"_l$(Par.l)"*arg)
generateFileName(Par::ParquetParams,arg::String = "") = _generateFileName(Par,"_p"*arg)

## ______________ State Variables shorthand ______________

"""Struct containing all memory used in a single ODE step """
struct MultiLoopWorkSpace{T,Buff,P <: PMFRGParams}
    State::StateType{T}
    Deriv::StateType{T}
    
    X::BubbleType{T} #full bubble
    XLeft::BubbleType{T} #left bubble
    XRight::BubbleType{T} #right bubble

    Y::BubbleType{T} #full bubble
    YLeft::BubbleType{T} #left bubble
    YRight::BubbleType{T} #right bubble

    Buffer::Buff #Buffer Arrays
    Par::P
end

function MultLoopWorkspace(Deriv::ArrayPartition,State::ArrayPartition,X,Y,Buffer,Par)
    error("Todo")
end

struct ParquetWorkspace{T,Buff,P <: PMFRGParams} <:PMFRGWorkspace
    State::StateType{T} #Stores the current step of iteration
    OldState::StateType{T} #Stores the last step of iteration

    Γ0::BareVertexType{T} #Stores the bare vertex
    I::VertexType{T} #Stores the irreducible vertex
    
    X::BubbleType{T} #Stores the actual bubble function X = B0 + BX

    B0::BubbleType{T} # Bubble involving bare vertex Γ_0 ∘ Γ
    BX::BubbleType{T} # Bubble involving other bubble X ∘ Γ

    Buffer::Buff #Buffer Arrays

    Par::P # Params
end

"""Constructs Workspace for parquet equations. Irreducible vertex can later be provided by optional argument but this is currently not supported """
function SetupParquet(Par::ParquetParams)
    @unpack NUnique,Npairs = Par.System
    IrrVertex = BareVertex_Freq(Par)

    PropsBuffers = [Matrix{double}(undef,NUnique,NUnique) for _ in 1:Threads.nthreads()] 
    VertexBuffers = [VertexBufferType(Npairs) for _ in 1:Threads.nthreads()]
    BubbleBuffers = [BubbleBufferType(Npairs) for _ in 1:Threads.nthreads()]
    Buffs = BufferTypeTwoLoop(PropsBuffers,VertexBuffers,BubbleBuffers) 
 
    Workspace = ParquetWorkspace(
        StateType(Par),
        StateType(Par),

        BareVertexType(Par),
        IrrVertex,
        
        BubbleType(Par),

        BubbleType(Par),
        BubbleType(Par),
        
        Buffs,
        Par
        )
    setToBareVertex!(Workspace.OldState.Γ,Par)
    setToBareVertex!(Workspace.State.Γ,Par)
    return Workspace
end
