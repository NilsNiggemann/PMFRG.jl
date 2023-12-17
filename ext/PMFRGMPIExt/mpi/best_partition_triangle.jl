module BestPartitionTriangle
include("./partition.jl")
import .Partition: partitions

struct TooManyRanksError <: Exception
    rank_tu::Int
    N::Int
end

function Base.show(io::IO, e::TooManyRanksError)
    println(io, "Too many ranks ($(e.rank_tu)) for N=$(e.N)")
end

"""Returns the t and u ranges for a given rank out of nranks,
   trying to produce the most balanced partition
   (where the weight is the number of sites of the given parity
   in the triangle where iu <= it).
"""
function _get_ranges_tu(N::Int, nranks_tu::Int, rank_tu::Int)
    nsites = get_number_of_sites_in_triangle(N)

    target_nsites_per_rank = nsites / nranks_tu

    nsites_before_target = round(target_nsites_per_rank * rank_tu)
    nsites_after_target = round(target_nsites_per_rank * (rank_tu + 1))

    t_start =
        argmin(x -> abs(get_number_of_sites_in_triangle(x) - nsites_before_target), 0:N) + 1
    t_stop = argmin(x -> abs(get_number_of_sites_in_triangle(x) - nsites_after_target), 1:N)

    if t_start > N
        throw(TooManyRanksError(rank_tu, N))
    end

    u_start, u_stop = 1, t_stop

    t_start:t_stop, u_start:u_stop

end


function _count_sites(itrange, iurange, isrange, parity)
    @assert iurange.start == 1
    @assert iurange.stop == itrange.stop
    function nsites_even_srange(itrange, isrange)
        (itrange.start + itrange.stop) * length(itrange) / 2 * length(isrange) / 2
    end
    if length(isrange) % 2 == 0
        nsites_even_srange(itrange, isrange)
    else
        # tu slice
        triangle_side = itrange.stop - itrange.start + 1
        # rectangle
        sir_a = div(triangle_side * itrange.start, 2, RoundFromZero)
        sir_b = div(triangle_side * itrange.start, 2, RoundToZero)

        if triangle_side % 2 == 1 && itrange.start % 2 == 1

            parity_of_corners = (itrange.start + 1) % 2
            sir_e, sir_o = [sir_a, sir_b][[parity_of_corners + 1, 2 - parity_of_corners]]
        else
            sir_e = sir_a
            sir_o = sir_b
        end
        # triangle
        # sit_e, sit_o: sites in triangle (t+u even/odd)
        sit_e, sit_o = _get_number_of_sites_eo(triangle_side - 1)

        # sip_e, sip_o: sites in traPeze (t+u even/odd)
        sip_e = sir_e + sit_e
        sip_o = sir_o + sit_o

        nsites_base = nsites_even_srange(itrange, isrange.start:(isrange.stop-1))
        #nsites_base + [sip_e,sip_o][parity + 1] # NOPE
        eo_parity = (isrange.stop + parity) % 2
        nsites_base + [sip_e, sip_o][1+eo_parity]
    end
end


function get_imbalance_from_ranges(N::Int, nranks::Int, all_ranges, parity::Int)
    min_nsites = typemax(Int64)
    max_nsites = 0
    for (irank, (isrange, itrange, iurange)) in enumerate(all_ranges)
        nsites = _count_sites(itrange, iurange, isrange, parity)
        min_nsites = (nsites < min_nsites) ? nsites : min_nsites
        max_nsites = (nsites > max_nsites) ? nsites : max_nsites
    end
    return (max_nsites - min_nsites) / max_nsites
end



"""Returns an estimate of the load imbalance among the ranks,
   in the range 0.0-1.0.
 """
function get_imbalance(N, nranks, get_ranges_func, parity)
    all_ranges = get_ranges_func(N, nranks, parity)
    get_imbalance_from_ranges(N, nranks, all_ranges, parity)
end


function _all_ns_x_ntu_factorizations(nranks)
    [(nranks_s, div(nranks, nranks_s)) for nranks_s = 1:nranks if (nranks % nranks_s) == 0]
end

"""
For given values of N, nranks and parity choices,
returns the ranges in s,t,u that correspond to the best balance partition,
for all the ranks.
"""
function get_all_ranges_stu(N::Int, nranks::Int, parity::Int)
    imbalance = 1.0

    all_ranges_best = Vector{Tuple{UnitRange{Int64},UnitRange{Int64},UnitRange{Int64}}}()

    for (nranks_s, nranks_tu) in _all_ns_x_ntu_factorizations(nranks)
        all_ranges = Vector{Tuple{UnitRange,UnitRange,UnitRange}}()
        for rank = 0:(nranks-1)

            rank_s = rank % nranks_s
            range_s = partitions(N, nranks_s)[1+rank_s]

            rank_tu = div(rank, nranks_s, RoundToZero)
            range_tu = _get_ranges_tu(N, nranks_tu, rank_tu)

            range_stu = (range_s, range_tu...)
            push!(all_ranges, range_stu)
        end
        candidate_imbalance = get_imbalance_from_ranges(N, nranks, all_ranges, parity)
        if candidate_imbalance < imbalance
            all_ranges_best = all_ranges
            imbalance = candidate_imbalance
        end
    end
    all_ranges_best
end

"""Returns the number of sites having iu<=it<=Ntu"""
function get_number_of_sites_in_triangle(Ntu)
    Int(Ntu * (Ntu + 1) / 2)
end

"""Returns the number of sites having iu<=it<=Ntu with a given parity."""
function _get_number_of_sites_eo(Ntu)
    total_elements = get_number_of_sites_in_triangle(Ntu)
    even = if (Ntu % 2 == 0)
        nhalf = Ntu / 2
        2 * nhalf * (nhalf + 1) / 2
    else
        nhalf = (Ntu - 1) / 2
        2 * nhalf * (nhalf + 1) / 2 + nhalf + 1
    end

    odd = total_elements - even
    even, odd
end
end
