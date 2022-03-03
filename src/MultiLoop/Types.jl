
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

function setToBareVertex!(Γ::VertexType,couplings::AbstractVector)
    Γ.a .= 0.
    Γ.b .= 0.
    for Rj in eachindex(couplings,axes(Γ.c,1))
        Γ.c[Rj,:,:,:] .= -couplings[Rj]
    end
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
