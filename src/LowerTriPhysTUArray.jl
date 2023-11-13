module LowerTriPhysTUArray
export Even, Odd, lowerTriPhysTUArray

@enum Parity Even Odd

"""A 4-dimensional array of sizes [Npairs,Ns,Ntu,Ntu]
   that contains only the elements for which:
   - iu <= it (or it <= iu)
   - is+it+iu has the given parity.
   Whether the array stores
   the elements for which iu <= it
   or the elements for which it <= iu
   depends on the implementation of the indexing functions.
"""
struct lowerTriPhysTUArray{T} <: AbstractArray{T,4}
    v::Vector{T} # The 3rd and 4rd dimensions are unrolled

    Npairs::Int64 # Int32 should be more than enough
    # s_extent and tu_extent very likely will be equal,
    # keeping this more general just in case it is needed later.
    s_extent::Int64 # Int32 should be more than enough
    tu_extent::Int64 # Int32 should be more than enough
    length::Int64
    parity::Int64 # Needs only 0 or 1 actually
end

function lowerTriPhysTUArray{T}(Npairs::Int,
                            s_extent::Int,
                            tu_extent::Int,
                            parity::Parity) where T

    eo_arr = [n_eo_elements_in_tri_half_stu(s_extent,tu_extent)...]
    parity_idx = if parity == Even 1 else 2 end

    length = eo_arr[parity_idx]*Npairs
    lowerTriPhysTUArray(Vector{T}(undef,length),
                    Npairs,
                    s_extent,
                    tu_extent,
                    length,
                    if parity == Even 0 else 1 end)
end

function Base.size(A::lowerTriPhysTUArray{T}) where T
    A.Npairs,A.s_extent,A.tu_extent,A.tu_extent
end

function Base.length(A::lowerTriPhysTUArray{T}) where T
    A.length
end

"""Function to convert at LowerTriPhysTUArray
   to a full array, for printing (and REPL testing)
"""
function full_repr(arr::lowerTriPhysTUArray{T}) where T
    full_repr = Array{Union{T,Missing},4}(undef,(arr.Npairs,
                                                 arr.s_extent,
                                                 arr.tu_extent,
                                                 arr.tu_extent))
    for idx in eachindex(IndexCartesian(),full_repr)
        rij,is,it,iu = Tuple(idx)
        try
            full_repr[rij,is,it,iu] = arr[rij,is,it,iu]
        catch e
            if e isa ParityError || e isa TriHalfError
                full_repr[rij,is,it,iu] = Missing()
            else
                throw(e)
            end
        end
    end
    full_repr
end

"""Method to display a LowerTriPhysTUArray in the REPL"""
function Base.show(io::IOContext{Base.TTY},
                   ::MIME{Symbol("text/plain")},
                   arr::lowerTriPhysTUArray{T}) where T
    unused = MIME{Symbol("text/plain")}()
    show(io,unused,full_repr(arr))
end

"""Generic show method (more testable)"""
function Base.show(io::IO,
                   arr::lowerTriPhysTUArray{T}) where T
    show(io,full_repr(arr))
end


"""Number of even and odd elements in the lower half
   (including the diagonal) of a square array of size Ntu x Ntu
"""
function n_eo_elements_in_lower_half_tu(Ntu::Int)
    nhalf = Int(floor(Ntu/2))
    total_elements = Int(Ntu*(Ntu+1)/2)

    a = nhalf*(nhalf+1)
    b = total_elements - a

    even, odd = if Ntu%2 == 0 a,b else b,a end
    even, odd
end

"""Number of even and odd elements in the lower half of the tu-plane
   (including the diagonal plane) of a Ns x Ntu x Ntu array
"""
function n_eo_elements_in_tri_half_stu(Ns::Int,Ntu::Int)
    all_lower_tu_half = Int(Ntu*(Ntu+1)/2)
    Ns_half = Int(floor(Ns/2))
    common = Ns_half*all_lower_tu_half
    if Ns%2 == 0
        e = o = common
        e,o
    else
        even_tu, odd_tu = n_eo_elements_in_lower_half_tu(Ntu)
        # The parity in the remaining s-odd slice is inverted
        even = common + odd_tu
        odd = common + even_tu
        even,odd
    end
end

# From https://docs.juliahub.com/Exceptions/IC1nl/0.1.0/manual/guide/
mutable struct TriHalfError <: Exception
    it::Int
    iu::Int
end
function Base.showerror(io::IO, e::TriHalfError)
    print(io,"TriangularHalfError:")
    print(io," Data not stored for index it > iu ($(e.it) > $(e.iu))")
end

"""Function used to compute the OFFSET in the [s,t,u] sector of a given element
   in a LowerTriPhysTUArray.
   It is assumed that iu <= it, otherwise an exception will be thrown.
"""
function _stu_eo_offset(s_extent::Int, is::Int, it::Int, iu::Int)
    if iu <= it
        a_offset = iu - 1
        b_offset = it - 1
    else
        throw(TriHalfError(it,iu))
    end
    s_offset = is - 1
    parity = (is+it+iu)%2

    e,o = n_eo_elements_in_tri_half_stu(s_extent,b_offset)

    eoarr = [e,o]
    return eoarr[parity+1] + Int(floor((a_offset*s_extent+s_offset)/2))

end

"""Function used to compute the OFFSET in the [Rij, s,t,u] sector of a given element
   in a LowerTriPhysTUArray.
"""
function _rstu_eo_offset(Npairs::Int, rij::Int,s_extent::Int, is::Int, it::Int, iu::Int)
    rij_offset = rij - 1
    _stu_eo_offset(s_extent,is,it,iu)*Npairs + rij_offset
end



"""
    _rstu_eo_idx(Npairs::Int, rij::Int,s_extent::Int, is::Int, it::Int, iu::Int)

Function used to compute the INDEX in the [Rij, s,t,u] sector of a given element
   in a LowerTriPhysTUArray.
"""
function _rstu_eo_idx(Npairs::Int, rij::Int,s_extent::Int, is::Int, it::Int, iu::Int)
    _rstu_eo_offset(Npairs, rij,s_extent, is, it, iu)+1
end


function get_rij(Npairs::Int,rstu_eo_idx::Int)
    rstu_eo_offset = rstu_eo_idx-1
    rij_offset = rstu_eo_offset % Npairs
    rij_offset + 1
end

function get_istu(Npairs::Int, s_extent::Int,  rstu_eo_idx::Int )
    1,1,1 # FIXME
end

mutable struct ParityError{A,I} <: Exception
    a::A
    i::I
end

function Base.showerror(io::IO, e::ParityError)
    parity = (e.i[1]+e.i[2]+e.i[3])%2
    print(io,"ParityError: data not stored for $(e.i), parity $(parity) ")
    println(io,"in array with parity $(e.a.parity):")
    Base.show(io,e.a)
end

function _check_stu_parity(A::lowerTriPhysTUArray{T}, is::Int,it::Int,iu::Int) where T
    if mod(is+it+iu, 2) != A.parity
        throw(ParityError(A,[is,it,iu]))
    end
end

function Base.setindex!(A::lowerTriPhysTUArray{T},
                        val::S,
                        Rij::Int64,
                        is::Int64,
                        it::Int64,
                        iu::Int64) where {T,S}
    _check_stu_parity(A,is,it,iu)
    A.v[_rstu_eo_idx(A.Npairs,Rij,A.s_extent,is,it,iu)] = T(val)
end

function Base.getindex(A::lowerTriPhysTUArray{T},
                       Rij::Int64,
                       is::Int64,
                       it::Int64,
                       iu::Int64) where T
    _check_stu_parity(A,is,it,iu)
    A.v[_rstu_eo_idx(A.Npairs,Rij,A.s_extent,is,it,iu)]
end
end
