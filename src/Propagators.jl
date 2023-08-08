θ(T) = 1 / sqrt(T)
Θ(ω,Λ,T) = 1/(1 + Λ^2/(T^2*ω^2))

function iG0(nw, T, Λ = 0.)
    w = get_w(nw)
    return θ(T) * Θ(w,Λ,T) / w
end
"""
Taking a Matsubara integer, gives fully dressed propagator
"""
function iG_(gamma::AbstractArray, x::Integer, T::Real, nw::Integer,Λ = 0.)
    w = get_w(nw)

    return 1/(√T *w + Λ^2/(T^1.5*w) + gamma_(gamma, x, nw))
end

"""
Taking a Matsubara integer, gives single-scale propagator without Katanin (self energy flow equation)
"""
function iS_(gamma::AbstractArray, x::Integer, T::Real, nw::Integer, Λ::Real = 0.)
    w = get_w(nw)
    # return -θ(T) * w / 2 * iG_(gamma, x, T, nw)^2
    return -(T^2*w^2 - 3*Λ^2)/(2. *T^2.5*w) * iG_(gamma, x, T, nw)^2
end

"""
Taking a Matsubara integer, gives full single-scale propagator with Katanin (vertex flow equation)
"""
function iSKat_(
    gamma::AbstractArray,
    Dgamma::AbstractArray,
    x::Integer,
    T::Real,
    nw::Integer,
    Λ::Real = 0.,
)
    w = get_w(nw)
    return -((T^2*w^2 - 3*Λ^2)/(2. *T^2.5*w) + gamma_(Dgamma, x, nw)) * iG_(gamma, x, T, nw)^2
end

"""given a Matsubara (!) integer, return the corresponding Matsubara frequency divided by temperature"""
function get_w(nw)
    return (2 * nw + 1) * pi
end

