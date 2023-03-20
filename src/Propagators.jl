θ(T) = 1/sqrt(T)
function iG0(nw,T)
    w = get_w(nw)
    θ(T)/(w)
end
"""
Taking a Matsubara integer, gives fully dressed propagator
"""
function iG_(gamma::AbstractArray, x::Integer, T::Real, nw::Integer)
    w = get_w(nw)
    return 1/(w/θ(T) + gamma_(gamma,x,nw))
end

"""
Taking a Matsubara integer, gives single-scale propagator without Katanin (self energy flow equation)
"""
function iS_(gamma::AbstractArray, x::Integer, T::Real, nw::Integer)
    w = get_w(nw)
    return -θ(T)*w/2 * iG_(gamma,x,T,nw)^2
end

"""
Taking a Matsubara integer, gives full single-scale propagator with Katanin (vertex flow equation)
"""
function iSKat_(gamma::AbstractArray, Dgamma::AbstractArray, x::Integer, T::Real, nw::Integer)
    w = get_w(nw)
    return -(θ(T)*w/2 + gamma_(Dgamma,x,nw) ) * iG_(gamma,x,T,nw)^2
end

"""given a Matsubara (!) integer, return the corresponding Matsubara frequency divided by temperature"""
function get_w(nw)
    return (2*nw+1)*pi
end

