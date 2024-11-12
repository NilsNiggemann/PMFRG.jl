"""Set two objects equal by recusively setting all their fields to be equal"""
function writeTo!(A::T, B::T) where {T<:VertexOrBubble}
    for f in fieldnames(T)
        a = getfield(A, f)
        b = getfield(B, f)
        writeTo!(a, b)
    end
    return A
end

function writeTo!(A::AbstractArray, B::AbstractArray)
    A .= B
    return A
end

"""Computes total difference between two objects of type T by summing up the distance of all elements"""
function squareDist(Γ1::T, Γ2::T) where {T<:VertexOrBubble}
    d = 0.0
    for f in fieldnames(T)
        a = getfield(Γ1, f)
        b = getfield(Γ2, f)
        d += squareDist(a, b)
    end
    return d
end

function squareDist(a::AbstractArray, b::AbstractArray)
    d = 0.0
    for i in eachindex(a, b)
        d += squareDist(a[i], b[i])
    end
    return d
end

@inline squareDist(a::Number, b::Number) = abs2(a - b)
@inline dist(a, b) = sqrt(squareDist(a, b))

@inline reldist(A, B) = dist(A, B) / max(norm(A), norm(B))

import SpinFRGLattices: squareNorm, norm

function squareNorm(A::T) where {T<:VertexOrBubble}
    s = 0.0
    for f in fieldnames(T)
        a = getfield(A, f)
        s += SpinFRGLattices.squareNorm(a)
    end
    return s
end

function norm(A)
    return sqrt(squareNorm(A))
end

function Threadsfill!(Tensor::AbstractArray, val)
    Threads.@threads for i in eachindex(Tensor)
        Tensor[i] = val
    end
end

"""fills an array with zeros"""
setZero!(a::AbstractArray{T,N}) where {T,N} = @tturbo a .= zero(T)

"""Recursively sets structure to zero"""
function setZero!(a::T) where {T<:VertexOrBubble}
    for f in fieldnames(T)
        setZero!(getfield(a, f))
    end
    return a
end
function absmax(x)
    MaxAndPos = findmax(x)
    MinAndPos = findmin(x)
    absvals = (abs(MaxAndPos[1]), abs(MinAndPos[1]))
    return (MaxAndPos, MinAndPos)[argmax(absvals)]
end
strd(x) = string(round(x, digits = 3))
