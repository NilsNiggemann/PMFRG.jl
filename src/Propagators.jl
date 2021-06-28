"""
Taking a Matsubara integer, gives fully dressed propagator
"""
function iG_(gamma,x::Integer, Lam,nw,Par)
    @unpack T = Par
    w = get_w(nw,T)
    return(w/(w^2+w*gamma_(gamma,x,nw,Par) + Lam^2))
end

"""
Taking a Matsubara integer, gives single-scale propagator without Katanin (self energy flow equation)
"""
function iS_(gamma,x::Integer, Lam,nw,Par)
    @unpack T= Par
    w = get_w(nw,T)
    return(-iG_(gamma,x,Lam,nw,Par)^2* 2*Lam/w )
end

"""
Taking a Matsubara integer, gives full single-scale propagator with Katanin (vertex flow equation)
"""
function iSKat_(gamma, Dgamma,x::Integer,Lam,nw,Par)
    @unpack T = Par
    w = get_w(nw,T)
    return(-iG_(gamma,x,Lam,nw,Par)^2*(2*Lam/w + gamma_(Dgamma,x,nw,Par)) )
end