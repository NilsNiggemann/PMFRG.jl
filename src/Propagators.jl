"""
Taking a Matsubara integer, gives fully dressed propagator
"""
function iG_(gamma,x::Integer, Lam,nw,Par)
    T = Par.NumericalParams.T
    w = get_w(nw,T)
    return(w/(w^2+w*gamma_(gamma,x,nw,Par) + Lam^2))
end

"""
Taking a Matsubara integer, gives single-scale propagator without Katanin (self energy flow equation)
"""
function iS_(gamma,x::Integer, Lam,nw,Par)
    T = Par.NumericalParams.T
    w = get_w(nw,T)
    return(-iG_(gamma,x,Lam,nw,Par)^2* 2*Lam/w )
end

"""
Taking a Matsubara integer, gives full single-scale propagator with Katanin (vertex flow equation)
"""
function iSKat_(gamma, Dgamma,x::Integer,Lam,nw,Par)
    T = Par.NumericalParams.T
    w = get_w(nw,T)
    return(-iG_(gamma,x,Lam,nw,Par)^2*(2*Lam/w + gamma_(Dgamma,x,nw,Par)) )
end

"""given a Matsubara (!) integer, return the corresponding Matsubara frequency"""
function get_w(nw,T)
    return (2*nw+1)*pi*T
end

