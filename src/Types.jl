"""Abstract Param struct to dispatch between different PMFRG methods (i.e. Two-Loop or Parquet)"""
abstract type PMFRGParams end

"""
Struct to hold all relevant quantities that are needed throughout the computation. 
    N::Int = 24
    Ngamma::Int = 100 #Number of gamma frequencies
    VDims::NTuple{4,Int} = (Npairs,N,N,N)
    accuracy::double = 1e-6
    Lam_min::double = 0.0
    Lam_max::double = 100.0
    usesymmetry::Bool = true
    MinimalOutput::Bool = false
    ex_freq::double = (2*N-1)*pi*T
    # lenIntw::Int = floor(Int,  (min(ex_freq,Lam_max)/pi/T-1)/2)
    lenIntw::Int = N
    lenIntw_acc::Int = Ngamma # more accurate for less demanding sums
    np_vec::Array{Int,1} = collect(0:N-1)
    np_vec_gamma::Array{Int,1} = collect(0:Ngamma-1)

"""
@proto struct NumericalParams{F} where {F <: AbstractFloat}
    T::double # Temperature
    N::Int # Number of positive frequencies
    Ngamma::Int  # Number of postivei gamma frequencies

    accuracy::F # convert type to float type
    Lam_min::F
    Lam_max::F
    
    lenIntw::Int 
    lenIntw_acc::Int  # more accurate for less demanding sums
    np_vec::Array{I,1}
    np_vec_gamma::Array{I,1}

    VDims::NTuple{4,I} # Vertex dimensions
    ex_freq::F
end

function NumericalParams(;
    T::F = 0.5, # Temperature
    N::Integer = 24,
    Ngamma::Integer = N, #Number of gamma frequencies
    accuracy::F = F(1e-6), # convert type to float type
    Lam_min::F = F(0.0),
    Lam_max::F = F(100.0),
    lenIntw::Int = N,
    lenIntw_acc::Int = Ngamma, # more accurate for less demanding sums
    np_vec::Array{I,1} = collect(0:N-1),
    np_vec_gamma::Array{I,1} = collect(0:Ngamma-1),
    kwargs...) where F <: AbstractFloat
    return NumericalParams(
        N,
        Ngamma ,
        accuracy,
        Lam_min,
        Lam_max ,
        usesymmetry ,
        MinimalOutput ,
        lenIntw = N,
        lenIntw_acc = Ngamma,
        np_vec,
        np_vec_gammanp_vec,
        (Npairs,N,N,N), #VDims,
        ex_freq = (2*N-1)*pi*T,
    )

abstract type AbstractOptions end 

@proto @with_kw struct OptionParams <: AbstractOptions
    usesymmetry::Bool = true
    MinimalOutput::Bool = false
end

"""Struct containing information about four-point vertices"""
struct VertexType{T}
    a::Array{T,4}
    b::Array{T,4}
    c::Array{T,4}
end


function VertexType(VDims::Tuple)
    return VertexType(
        zeros(VDims), # Gamma_a
        zeros(VDims), # Gamma_b
        zeros(VDims) # Gamma_c
    )
end

function setToBareVertex!(Γc::AbstractArray{T,4},couplings::AbstractVector) where T
    for Rj in eachindex(couplings,axes(Γc,1))
        Γc[Rj,:,:,:] .= -couplings[Rj]
    end
    return Γc
end

function setToBareVertex!(Γ::VertexType,couplings::AbstractVector)
    Γ.a .= 0.
    Γ.b .= 0.
    setToBareVertex!(Γ.c,couplings)
    return Γ
end

setToBareVertex!(Γ::VertexType,Par) = setToBareVertex!(Γ,Par.System.couplings)

"""Allocates a frequency-dependent vertex type and initializes it with the bare vertex"""
BareVertex_Freq(Par) = setToBareVertex!(VertexType(Par.VDims),Par.System.couplings)

"""Stores all information about the bare vertex (without frequencies)"""
struct BareVertexType{T}
    c::Vector{T}
end

BareVertexType(Par::Params) = BareVertexType(-Par.couplings)

"""Struct containing information about the (physical) ODE State, i.e. vertices"""
struct StateType{T}
    f_int::Array{T,1}
    γ::Array{T,2}
    Γ::VertexType{T}
end

function StateType(NUnique::Int,Ngamma::Int,VDims::Tuple)
    return StateType(
        zeros(NUnique), # fint
        zeros(NUnique,Ngamma), # gamma
        VertexType(VDims)
    )
end

StateType(Par::Params) = StateType(Par.NUnique,Par.Ngamma,Par.VDims) 

"""Struct storing information about Bubble functions, i.e. rhs derivatives of vertex flow equations."""
struct BubbleType{T}
    a::Array{T,4} #"a" type bubble
    b::Array{T,4}
    c::Array{T,4}

    Ta::Array{T,4} #"a-Tilde" type bubble
    Tb::Array{T,4}
    Tc::Array{T,4}
    Td::Array{T,4}
end

function BubbleType(VDims::Tuple)
    return BubbleType(
        zeros(VDims),
        zeros(VDims),
        zeros(VDims),

        zeros(VDims),
        zeros(VDims),
        zeros(VDims),
        zeros(VDims)
    )
end

BubbleType(Par::Params) = BubbleType(Par.VDims) 
    
struct VertexBufferType{T}
	Va12::Vector{T}
	Vb12::Vector{T}
	Vc12::Vector{T}

	Va34::Vector{T}
	Vb34::Vector{T}
	Vc34::Vector{T}
	
	Vc21::Vector{T}
	Vc43::Vector{T}
end
VertexBufferType(Npairs) = VertexBufferType((zeros(Npairs) for _ in 1:8)...)

struct BufferType{PropsBuff,VertexBuff}
    Props::Vector{PropsBuff}
    Vertex::Vector{VertexBuff}
end
