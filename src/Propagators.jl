"""
Taking a Matsubara integer, gives fully dressed propagator
"""
function iG(gamma,x::Integer, Lam,nw,Par)
    @unpack T = Par
    w = get_w(nw,T)
    return(w/(w^2+w*gamma_(gamma,x,nw,Par) + Lam^2))
end

"""
Taking a Matsubara integer, gives single-scale propagator without Katanin (self energy flow equation)
"""
function iS(gamma,x::Integer, Lam,nw,Par)
    @unpack T= Par
    w = get_w(nw,T)
    return(-iG(gamma,x,Lam,nw,Par)^2* 2*Lam/w )
end

"""
Taking a Matsubara integer, gives full single-scale propagator with Katanin (vertex flow equation)
"""
function iSKat(gamma, Dgamma,x::Integer,Lam,nw,Par)
    @unpack T = Par
    w = get_w(nw,T)
    return(-iG(gamma,x,Lam,nw,Par)^2*(2*Lam/w + gamma_(Dgamma,x,nw,Par)) )
end