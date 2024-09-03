"""Given a (Matsubara) integer, returns a corresponding signed integer within the considered range"""
function get_sign_iw(nw::Integer, N::Integer)
    s = sign(nw)
    nw_bounds = min(abs(nw), N - 1)
    return s * nw_bounds + 1
end

"""Returns self energy evaluated at the appropriate frequency"""
function gamma_(gamma::AbstractArray, x::Integer, nw::Integer)
    Ngamma = size(gamma, 2)
    s = 1
    if nw < 0
        nw = -nw - 1
        s = -1
    end
    iw = get_sign_iw(nw, Ngamma)
    return s * gamma[x, iw]
end

@inline function convertFreqArgs(ns, nt, nu, Nw)
    # @assert (ns+nt+nu) %2 != 0 "trying to convert wrong freqs $ns + $nt +  $nu = $(ns+nt+nu)"
    swapsites = nt * nu < 0
    ns, nt, nu = abs.((ns, nt, nu))
    ns = min(ns, Nw - 1 - (ns + Nw - 1) % 2)
    nt = min(nt, Nw - 1 - (nt + Nw - 1) % 2)
    nu = min(nu, Nw - 1 - (nu + Nw - 1) % 2)

    return ns, nt, nu, swapsites
end

function shrink_OLD(ns::Int, nt::Int, nu::Int, Nw::Int)
    # @assert (ns+nt+nu) %2 != 0 "$ns + $nt +  $nu = $(ns+nt+nu)"
    max_abs_stu = max(abs(ns), abs(nt), abs(nu))
    # if outside range, shrink (ns,nt,nu) and put back initial parity, but stay in same octant
    if max_abs_stu > Nw - 1
        ns_p = Int(round(ns * (Nw - 1.0) / max_abs_stu))
        nt_p = Int(round(nt * (Nw - 1.0) / max_abs_stu))
        nu_p = Int(round(nu * (Nw - 1.0) / max_abs_stu))
        ns_p = sign(ns) * abs(ns_p - (ns_p + ns) % 2)
        nt_p = sign(nt) * abs(nt_p - (nt_p + nt) % 2)
        nu_p = sign(nu) * abs(nu_p - (nu_p + nu) % 2)
        return ns_p, nt_p, nu_p
    else
        return 1 * ns, 1 * nt, 1 * nu
    end
end
"""
Returns value of vertex, swaps sites i <-> j by reading from inverted pairs when necessary. For performance, it is advised to directly use the Arrays like Va after pre-checking the frequency structure 
"""
@inline function V_(
    Vertex::AbstractArray,
    Rj::Integer,
    ns::Integer,
    nt::Integer,
    nu::Integer,
    Rji::Integer,
    N::Integer,
)
    # @assert (ns+nt+nu) %2 != 0 "$ns + $nt +  $nu = $(ns+nt+nu)"
    ns, nt, nu, swapsites = convertFreqArgs(ns, nt, nu, N)
    Rj = ifelse(swapsites, Rji, Rj)
    return @inbounds Vertex[Rj, ns+1, nt+1, nu+1]
end

@inline function bufferV_!(
    Cache,
    Vertex::AbstractArray,
    ns::Integer,
    nt::Integer,
    nu::Integer,
    invpairs::AbstractArray,
    N,
)
    ns, nt, nu, swapsites = convertFreqArgs(ns, nt, nu, N)
    # @assert (ns+nt+nu) %2 != 0 "$ns + $nt +  $nu = $(ns+nt+nu)"

    is, it, iu = ns + 1, nt + 1, nu + 1
    @inbounds begin
        if swapsites
            @turbo unroll = 1 inline = true for R in eachindex(Cache, invpairs)
                Cache[R] = Vertex[invpairs[R], is, it, iu]
            end
        else
            @turbo unroll = 1 inline = true for R in eachindex(Cache, invpairs)
                Cache[R] = Vertex[R, is, it, iu]
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
