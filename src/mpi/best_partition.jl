module BestPartition
include("./partition.jl")
using .Partition: partitions
include("./decompose.jl")
using .Decompose: decompose

function imbalance(rangess::Vector{Vector{UnitRange{Int64}}})
    min_volume = prod(sum(length(range) for range in ranges) for ranges in rangess)
    max_volume = 0

    nranges = [length(ranges) for ranges in rangess]
    for idxs in Base.Iterators.product([1:n for n in nranges]...)
        vol = prod(length(ranges[idx]) for (ranges,idx) in zip(rangess,idxs))
        min_volume = min(min_volume,vol)
        max_volume = max(max_volume,vol)
    end

    max_volume - min_volume
end


function find_best_ND_partition(dimensions,ranks)
    nfactors = length(dimensions)
    rank_decompositions = decompose(ranks,nfactors)

    function get_imbalance(rank_decomposition)
        ranges = [ partitions(l,rd) for (l,rd) in zip(dimensions,rank_decomposition) ]
        imbalance(ranges)
    end
    best_decomposition = argmin(get_imbalance,rank_decompositions)
    #best_imbalance = get_imbalance(best_decomposition)
    best_decomposition

end

"""
Given a N-dimensional array, find the best decomposition of it
and then return the ranges on each dimension
that correspond to the current rank.
rank lies in the range [0,nranks)
"""
function get_ranges(dimensions,nranks,rank)
        decomposition = find_best_ND_partition(dimensions,nranks)

        idxs = collect(Base.Iterators.product([1:npart for npart in decomposition]...))[rank+1]

        iranges = [ partitions(L,npart)[idx]
                    for (L,npart,idx) in
                        zip(dimensions,decomposition,idxs)]

        iranges
end
end
