
@enum Parity Even Odd

"""A 4-dimensional array that is symmetric in the 2 last dimensions."""
struct SymmPhysTUArray{T} <: AbstractArray{T,4}
    v::Vector{T} # The 3rd and 4rd dimensions are unrolled

    Npairs::Int64 # Int32 should be more than enough
    # s_extent and tu_extent very likely will be equal,
    # keeping this more general just in case it is needed later.
    s_extent::Int64 # Int32 should be more than enough
    tu_extent::Int64 # Int32 should be more than enough
    length::Int64
    parity::Int64 # Needs only 0 or 1 actually
end



function SymmPhysTUArray{T}(Npairs::Int,
                            s_extent::Int,
                            tu_extent::Int,
                            parity::Parity) where T
    length = Npairs*Int(s_extent*tu_extent*(tu_extent+1) / 2 / 2) # ~+- 1 # FIXME
    SymmPhysTUArray(Vector{T}(undef,length),
                    Npairs,
                    s_extent,
                    tu_extent,
                    length,
                    if parity == Even 0 else 1 end)
end

function Base.size(A::SymmPhysTUArray{T}) where T
    A.Npairs,A.s_extent,A.tu_extent,A.tu_extent
end

function Base.length(A::SymmPhysTUArray{T}) where T
    A.length
end


function _tu_offset(it::Int, iu::Int)
    if it > iu
        a_offset = iu - 1
        b_offset = it - 1
    else
        a_offset = it - 1
        b_offset = iu - 1
    end
    return a_offset + Int((b_offset+1)*b_offset/2)
end

function n_eo_elements_in_lower_half(N::Int)
    nhalf = Int(floor(N/2))
    total_elements = Int(N*(N+1)/2)

    a = nhalf*(nhalf+1)
    b = total_elements - a

    even, odd = if N%2 == 0 a,b else b,a end
    even, odd
end

function _tu_eo_offset(it::Int, iu::Int)
    if it < iu
        a_offset = it - 1
        b_offset = iu - 1
    else
        a_offset = iu - 1
        b_offset = it - 1
    end
    parity = (iu + it)%2
    e,o = n_eo_elements_in_lower_half(b_offset)
    eoarr = [e,o]

    return eoarr[parity+1] + Int(floor(a_offset/2))
end

function _stu_offset(s_extent::Int, is::Int, it::Int, iu::Int)
    s_offset = is-1
    _tu_offset( it, iu)*s_extent+s_offset
end

function _stu_eo_offset(s_extent::Int, is::Int, it::Int, iu::Int)
    Int(floor(_stu_offset(s_extent,
                          is,
                          it,
                          iu) / 2))
end

function _rstu_eo_offset(Npairs::Int, rij::Int,s_extent::Int, is::Int, it::Int, iu::Int)
    rij_offset = rij - 1
    _stu_eo_offset(s_extent,is,it,iu)*Npairs + rij_offset
end


function _rstu_eo_idx(Npairs::Int, rij::Int,s_extent::Int, is::Int, it::Int, iu::Int)
    _rstu_eo_offset(Npairs, rij,s_extent, is, it, iu)+1
end

#function site_parity(is::Int,it::Int,iu::Int)
#    mod(is+it+iu, 2)
#end


function _check_stu_parity(A::SymmPhysTUArray{T}, is::Int,it::Int,iu::Int) where T
    if mod(is+it+iu, 2) != A.parity
        throw(BoundsError(A,[is,it,iu]))
    end
end


function Base.setindex!(A::SymmPhysTUArray{T},
                        val::S,
                        Rij::Int64,
                        is::Int64,
                        it::Int64,
                        iu::Int64) where {T, S<:T}
    _check_stu_parity(A,is,it,iu)
    A.v[_rstu_eo_idx(A.Npairs,Rij,A.s_extent,is,it,iu)] = val
end

function Base.getindex(A::SymmPhysTUArray{T},
                       Rij::Int64,
                       is::Int64,
                       it::Int64,
                       iu::Int64) where T
    _check_stu_parity(A,is,it,iu)
    A.v[_rstu_eo_idx(A.Npairs,Rij,A.s_extent,is,it,iu)]
end
