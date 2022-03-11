
"""Set two objects equal by recusively setting all their fields to be equal"""
function writeTo!(A::T,B::T) where T
    for f in fieldnames(T)
        a = getfield(A,f)
        b = getfield(B,f)
        writeTo!(a,b)
    end
    return
end

function writeTo!(A::AbstractArray,B::AbstractArray)
    A .= B
    return
end

"""Computes total difference between two objects of type T by summing up the distance of all elements"""
function dist(Γ1::T,Γ2::T) where T
    d = 0.
    for f in fieldnames(VertexType)
        a = getfield(Γ1,f)
        b = getfield(Γ2,f)
        d += dist(a,b)
    end
    return d
end

function dist(a::AbstractArray, b::AbstractArray)
    d = 0.
    for i in eachindex(a,b)
        d += abs(a[i]-b[i])
    end
    return d
end

function Threadsfill!(Tensor::AbstractArray,val)
    Threads.@threads for i in eachindex(Tensor)
        Tensor[i] = val
    end
end

"""fills an array with zeros"""
setZero!(a::AbstractArray{T,N}) where {T,N} = fill!(a,zero(T))

function setZero!(PartArr::ArrayPartition)
    for arr in PartArr.x
        fill!(arr,0.)
    end
end


"""Recursively sets structure to zero"""
function setZero!(a::T) where T
    for f in fieldnames(T)
        setZero!( getfield(a,f))
    end
    return a
end

