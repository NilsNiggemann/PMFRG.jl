struct MultiLoop
    l::Int
end

struct Parquet end

struct MultiLoopParams{F,G<:Geometry} <: PMFRGParams
    System::G
    NumericalParams::NumericalParams{F}
    Options::OptionParams
    l::Int
end

struct ParquetOptions{F<:AbstractFloat} <: AbstractOptions
    BSE_iters::Int
    SDE_iters::Int
    SDE_tolerance::F
    BSE_epsilon::F
    SDE_epsilon::F
    BSE_vel::F
    SDE_vel::F
    usesymmetry::Bool
    MinimalOutput::Bool
end

function ParquetOptions(;
    BSE_iters::Int = 40,
    SDE_iters::Int = 1000,
    SDE_tolerance::AbstractFloat = 1e-9,
    BSE_epsilon::AbstractFloat = 1.0,
    SDE_epsilon::AbstractFloat = 1.0,
    BSE_vel::AbstractFloat = 0.0,
    SDE_vel::AbstractFloat = 0.0,
    usesymmetry::Bool = true,
    MinimalOutput::Bool = false,
    kwargs...,
)
    return ParquetOptions(
        BSE_iters,
        SDE_iters,
        SDE_tolerance,
        BSE_epsilon,
        SDE_epsilon,
        BSE_vel,
        SDE_vel,
        usesymmetry,
        MinimalOutput,
    )
end

struct ParquetParams{F,G<:Geometry} <: PMFRGParams
    System::G
    NumericalParams::NumericalParams{F}
    Options::ParquetOptions{F}
end

Params(System::Geometry, O::MultiLoop; kwargs...) =
    MultiLoopParams(System, NumericalParams(; kwargs...), OptionParams(; kwargs...), O.l)

"""Params constructor for parquet calculation. For all available params, see also PMFRG.ParquetOptions and PMFRG.NumericalParams."""
function Params(System::Geometry, O::Parquet; eps = nothing, vel = 0.0, kwargs...)
    eps === nothing && (eps = getEpsilon(kwargs[:T]))
    PO = ParquetOptions(;
        BSE_epsilon = eps,
        SDE_epsilon = eps,
        BSE_vel = vel,
        SDE_vel = vel,
        kwargs...,
    )
    ParquetParams(System, NumericalParams(; kwargs...), PO)
end

function getEpsilon(T)
    min(T / 2, 1.0)
end

getLoopOrder(P::MultiLoopParams) = P.l
getPMFRGMethod(P::MultiLoopParams) = MultiLoop(getLoopOrder(P))
getPMFRGMethod(P::ParquetParams) = Parquet()

## ______________ State Variables shorthand ______________

"""Struct containing all memory used in a single ODE step """
struct MultiLoopWorkSpace{T,Buff,P<:PMFRGParams}
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

function MultLoopWorkspace(Deriv::AbstractVector, State::AbstractVector, X, Y, Buffer, Par)
    error("Todo")
end

struct ParquetWorkspace{T,Buff,P<:PMFRGParams} <: PMFRGWorkspace
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
    (; NUnique, Npairs) = Par.System
    IrrVertex = BareVertex_Freq(Par)

    floattype = _getFloatType(Par) #get type of float, i.e. Float64
    PropsBuffers = getChannel([
        MMatrix{NUnique,NUnique,floattype,NUnique * NUnique}(undef) for
        _ = 1:Threads.nthreads()
    ])
    VertexBuffers =
        getChannel([VertexBufferType(floattype, Npairs) for _ = 1:Threads.nthreads()])
    BubbleBuffers =
        getChannel([BubbleBufferType(floattype, Npairs) for _ = 1:Threads.nthreads()])
    Buffs = BufferTypeTwoLoop(PropsBuffers, VertexBuffers, BubbleBuffers)

    Workspace = ParquetWorkspace(
        StateType(Par),
        StateType(Par),
        BareVertexType(Par),
        IrrVertex,
        BubbleType(Par),
        BubbleType(Par),
        BubbleType(Par),
        Buffs,
        Par,
    )
    setToBareVertex!(Workspace.OldState.Γ, Par)
    setToBareVertex!(Workspace.State.Γ, Par)
    return Workspace
end
