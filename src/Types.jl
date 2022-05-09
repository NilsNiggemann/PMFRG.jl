const double = Float64

"""
Abstract Param struct to dispatch between different PMFRG methods (i.e. Two-Loop or Parquet)
Assumed to have at least the fields:
System::Geometry
NumericalParams::NumericalParams
Options
"""
abstract type PMFRGParams end
abstract type PMFRGWorkspace end

"""
Struct to hold all relevant quantities that are needed throughout the computation. 

    N::Int = 24
    Ngamma::Int = 100 #Number of gamma frequencies
    VDims::NTuple{4,Int} = (Npairs,N,N,N)
    accuracy::double = 1e-6
    Lam_min::double = 0.0
    Lam_max::double = 100.0 + 500*T,
    usesymmetry::Bool = true
    MinimalOutput::Bool = false
    ex_freq::double = (2*N-1)*pi*T
    lenIntw::Int = N
    lenIntw_acc::Int = 2*maximum((N,Ngamma,lenIntw)) # more accurate for less demanding sums
    np_vec::Array{Int,1} = collect(0:N-1)
    np_vec_gamma::Array{Int,1} = collect(0:Ngamma-1)

"""
struct NumericalParams{F <: AbstractFloat}
    T::F # Temperature
    N::Int # Number of positive frequencies
    Ngamma::Int  # Number of positive gamma frequencies

    accuracy::F # 
    Lam_min::F
    Lam_max::F
    
    lenIntw::Int 
    lenIntw_acc::Int  # more accurate for less demanding sums
    np_vec::Array{Int,1}
    np_vec_gamma::Array{Int,1}

    ex_freq::F
end

function NumericalParams(;
    T::AbstractFloat = 0.5, # Temperature
    N::Integer = 24,
    Ngamma::Integer = N, #Number of gamma frequencies
    accuracy::AbstractFloat = 1e-6, # convert type to float type
    Lam_min::AbstractFloat = 0.0,
    Lam_max::AbstractFloat = 100.0 + 500*T,
    lenIntw::Int = N,
    lenIntw_acc::Int = 2*maximum((N,Ngamma,lenIntw)), # more accurate for less demanding sums
    np_vec::Array{Int,1} = collect(0:N-1),
    np_vec_gamma::Array{Int,1} = collect(0:Ngamma-1),
    ex_freq = (2*N-1)*pi*T,
    kwargs...)
    return NumericalParams(
        T,
        N,
        Ngamma,
        accuracy,
        Lam_min,
        Lam_max ,
        lenIntw,
        lenIntw_acc,
        np_vec,
        np_vec_gamma,
        ex_freq,
    )
end
abstract type AbstractOptions end 
getPMFRGMethod(::Val{1}) = OneLoop()
getPMFRGMethod(::Val{2}) = TwoLoop()
getPMFRGMethod(::Val{n}) where n = MultiLoop(n)
getPMFRGMethod(n::Int) = getPMFRGMethod(Val(n))

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
BareVertex_Freq(Par) = setToBareVertex!(VertexType(getVDims(Par)),Par.System.couplings)

"""Stores all information about the bare vertex (without frequencies)"""
struct BareVertexType{T}
    c::Vector{T}
end

BareVertexType(Par::PMFRGParams) = BareVertexType(-Par.System.couplings)

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

StateType(Par::PMFRGParams) = StateType(Par.System.NUnique,Par.NumericalParams.Ngamma,getVDims(Par)) 

StateType(f_int,γ,Γa,Γb,Γc) = StateType(f_int,γ,VertexType(Γa,Γb,Γc)) 

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

const VertexOrBubble = Union{StateType,BubbleType,VertexType}

getVDims(Par::PMFRGParams) = (Par.System.Npairs,Par.NumericalParams.N,Par.NumericalParams.N,Par.NumericalParams.N)

BubbleType(Par::PMFRGParams) = BubbleType(getVDims(Par)) 
VertexType(Par::PMFRGParams) = VertexType(getVDims(Par)) 

function constructBubbleFromVertex!(B::BubbleType,Γ::VertexType)
    B.a .= Γ.a
    B.b .= Γ.b
    B.c .= Γ.c
    for Rij in axes(Γ.a,1), s in axes(Γ.a,2), t in axes(Γ.a,3), u in axes(Γ.a,4) 
        B.Ta[Rij,s,t,u] = -Γ.a[Rij,t,s,u] 
        B.Tb[Rij,s,t,u] = -Γ.c[Rij,t,s,u] 
        B.Tc[Rij,s,t,u] = -Γ.b[Rij,t,s,u]
        B.Td[Rij,s,t,u] = Γ.c[Rij,t,u,s]
    end
    return B
end
constructBubbleFromVertex(Γ::VertexType) = constructBubbleFromVertex!(BubbleType(size(Γ.a)),Γ) 

"""Struct containing the observables that are saved at every step"""
struct Observables{T}
    Chi::Vector{T}
    gamma::Matrix{T}
    f_int::Vector{T}
    MaxVa::Vector{T}
    MaxVb::Vector{T}
    MaxVc::Vector{T}
end


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
##
import RecursiveArrayTools:ArrayPartition

ArrayPartition(x::StateType) = ArrayPartition(x.f_int,x.γ,x.Γ.a,x.Γ.b,x.Γ.c) 