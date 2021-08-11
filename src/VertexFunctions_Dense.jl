##
"""Redefine sign function to include 0 as positive"""
function sgn(x::Number)
    if x >= 0
        return 1
    end
    return -1
end

"""Given a (Matsubara) integer, returns a corresponding signed integer within the considered range"""
function get_sign_iw(nw::Integer,N::Integer)
    s = sgn(nw)
    nw = abs(nw)
    if nw < N
        return s*(nw +1) # Need to add 1 because Julia arrays do not start at 0
    end
    return s *N
end

"""Given a (Matsubara) integer, returns a corresponding signed integer within the considered range"""
function get_sign_iw(nw::Integer,N::Integer)
    s = sgn(nw)
    nw = abs(nw)
    if nw < N
        return s*(nw +1) # Need to add 1 because Julia arrays do not start at 0
    end
    return s *N
end

"""Given the last two frequencies of a vertex return their absolute values together with a bool whether the pair ij has to be inverted"""
function flip_tu(it,iu)
    if it*iu>=0
        return (false, abs(it) , abs(iu))
    end
    return (true, abs(it) , abs(iu))
end

"""Returns self energy evaluated at the appropriate frequency"""
function gamma_(gamma::AbstractArray, x::Integer, nw::Integer,Par::Params)
    @unpack Ngamma = Par
    s = 1
    if nw<0
        nw = -nw -1
        s = -1
    end
    iw = get_sign_iw(nw,Ngamma)
    return s*gamma[x,iw]
end

@inline function convertFreqArgs(is,it,iu,N)
    swapsites = ifelse(it*iu<0,true, false)
    is,it,iu = abs.((is,it,iu))
    if (is+it+iu-3) %2 == 0
        is == N && return is-1,it,iu,swapsites
        it == N && return is,it-1,iu,swapsites
        iu == N && return is,it,iu-1,swapsites
    end
    @assert (is+it+iu-3) %2 != 0 "$is + $it +  $iu = $(is+it+iu)"
    return is,it,iu,swapsites
end

function shrink_OLD(ns::Int,nt::Int,nu::Int, Nw::Int)
    @assert (ns+nt+nu) %2 != 0 "$ns + $nt +  $nu = $(ns+nt+nu)"
    max_abs_stu = max(abs(ns),abs(nt),abs(nu))
    # if outside range, shrink (ns,nt,nu) and put back initial parity, but stay in same octant
    if max_abs_stu > Nw-1
        ns_p = Int(round(ns*(Nw-1.)/max_abs_stu))
        nt_p = Int(round(nt*(Nw-1.)/max_abs_stu))
        nu_p = Int(round(nu*(Nw-1.)/max_abs_stu))
        ns_p = sign(ns)*abs(ns_p - (ns_p+ns)%2)
        nt_p = sign(nt)*abs(nt_p - (nt_p+nt)%2)
        nu_p = sign(nu)*abs(nu_p - (nu_p+nu)%2)
        return ns_p,nt_p,nu_p
    else
        return 1*ns,1*nt,1*nu
    end
end
"""
Returns value of vertex, swaps sites i <-> j by reading from inverted pairs when necessary. For performance, it is advised to directly use the Arrays like Va after pre-checking the frequency structure 
"""
@inline function V_(Vertex::AbstractArray, Rj::Integer, is::Integer,it::Integer,iu::Integer,Rji::Integer,N::Integer)
    is,it,iu,swapsites = convertFreqArgs(is,it,iu,N)
    Rj = ifelse(swapsites,Rji,Rj)
    @inbounds Vertex[Rj,abs(is),abs(it),abs(iu)]
end

@inline function bufferV_!(Cache, Vertex::AbstractArray, is::Integer,it::Integer,iu::Integer,invpairs::AbstractArray,N)
    is,it,iu,swapsites = convertFreqArgs(is,it,iu,N)
    @inbounds begin 
        if swapsites<0
            @turbo unroll = 1 inline = true for R in eachindex(Cache,invpairs)
                Cache[R] = Vertex[invpairs[R],abs(is),abs(it),abs(iu)]
            end
        else
            @turbo unroll = 1 inline = true for R in eachindex(Cache,invpairs)
                Cache[R] = Vertex[R,abs(is),abs(it),abs(iu)]
            end
        end
    end
end

# function V_(Vertex::AbstractArray, Rj::Integer, is::Integer,it::Integer,iu::Integer,Rji::Integer)
#     if it*iu<0
#         return Vertex[Rji,abs(is),abs(it),abs(iu)]
#     end
#     return Vertex[Rj,abs(is),abs(it),abs(iu)]
# end