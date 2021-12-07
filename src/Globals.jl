const double = Float64
abstract type PMFRGMethod end
struct OneLoop <: PMFRGMethod end
"""
Struct to hold all relevant quantities that are needed throughout the computation. They can be unpacked using the Parameters package. 
Julia-Performance-tip: Make sure the type definitions are as accurate as possible to prevent run-time allocations. If used correctly, the values are only passed by reference.
"""
@with_kw struct Params{SType}
    # ______________ Model ______________

    T::double = 0.5 # Temperature

    # ______________ Geometry ______________
    System::Geometry
    Name::String = System.Name
    Npairs::Int = System.Npairs
    Nsum::Vector{Int} = System.Nsum
    NUnique::Int = System.NUnique
    couplings::Vector{double} = System.couplings
    invpairs::Vector{Int} = System.invpairs
    siteSum::Matrix{sumElements} = System.siteSum
    S::SType = StructArray(siteSum)
    PairTypes::Vector{sitePair} = System.PairTypes
    OnsitePairs::Vector{Int}= System.OnsitePairs
    
    # ______________ parameters for numerics ______________
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

end

"""Struct containing the observables that are saved at every step"""
struct Observables
    Chi::Vector{double}
    gamma::Matrix{double}
    f_int::Vector{double}
    MaxVa::Vector{double}
    MaxVb::Vector{double}
    MaxVc::Vector{double}
end

## ______________ State Variables shorthand ______________

"""
Convenience struct containing references to arrays for vertices and their derivative. This does not allocate any memory and is performant. 
"""
struct Workspace_Struct
    f_int::Array{double,1}
    gamma::Array{double,2}
    Va::Array{double,4}
    Vb::Array{double,4}
    Vc::Array{double,4}
    
    Df_int::Array{double,1}
    Dgamma::Array{double,2}
    DVa::Array{double,4}
    DVb::Array{double,4}
    DVc::Array{double,4}
    
    Xa::Array{double,4}
    Xb::Array{double,4}
    Xc::Array{double,4}
    XTa::Array{double,4}
    XTb::Array{double,4}
    XTc::Array{double,4}
    XTd::Array{double,4}
end

function Threadsfill!(Tensor,val)
    Threads.@threads for i in eachindex(Tensor)
        Tensor[i] = val
    end
end

function setZero!(PartArr)
    for arr in PartArr.x
        fill!(arr,0.)
    end
end

function Workspace_Struct(Deriv,State,X,XTilde)
    setZero!(Deriv)
    setZero!(X)
    setZero!(XTilde)
    return Workspace_Struct(State.x...,Deriv.x...,X.x...,XTilde.x...)
end
# given a Matsubara (!) integer, return the corresponding Matsubara frequency
function get_w(nw,T)
    return (2*nw+1)*pi*T
end

function absmax(x)
    MaxAndPos = findmax(x)
    MinAndPos = findmin(x)
    absvals = (abs(MaxAndPos[1]),abs(MinAndPos[1]))
    return (MaxAndPos,MinAndPos)[argmax(absvals)]
end
strd(x) = string(round(x,digits=3))