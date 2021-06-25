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

"""Given a (Matsubara) integer, returns a pair containing its sign and a corresponding positive integer within the considered range"""
function get_pos_iw(nw::Integer,N::Integer)
    s = sgn(nw)
    nw = abs(nw)
    if nw < N
        return s,nw+1
    end
    return s,N
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


"""
Returns value of vertex, swaps sites i <-> j by reading from inverted pairs when necessary. For performance, it is advised to directly use the Arrays like Va after pre-checking the frequency structure 
"""
function V_(Vertex::AbstractArray, Rj::Integer, is::Integer,it::Integer,iu::Integer,Rji::Integer)
    if it*iu<0
        return Vertex[Rji,abs(is),abs(it),abs(iu)]
    end
    return Vertex[Rj,abs(is),abs(it),abs(iu)]
end