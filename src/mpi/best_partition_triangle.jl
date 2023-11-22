module BestPartitionTriangle
include("./partition.jl")
import .Partition: partitions


"""Returns the t and u ranges for a given rank out of nranks,
   trying to produce the most balanced partition
   (where the weight is the number of sites of the given parity
   in the triangle where iu <= it).
"""
function _get_ranges_tu(N,nranks_tu, rank_tu)
    nsites = get_number_of_sites(N)
    nsites_per_rank = nsites/nranks_tu
    nsites_before_target = round(nsites_per_rank*rank_tu)
    nsites_after_target = round(nsites_per_rank*(rank_tu+1))

    start = argmin(x -> abs(get_number_of_sites(x) - nsites_before_target),0:N)
    end_ = argmin(x -> abs(get_number_of_sites(x) - nsites_after_target),1:N)

    start+1:end_ , 1:end_
end

function get_imbalance_from_ranges(N,nranks,all_ranges,parity)
    min_nsites = typemax(Int64)
    max_nsites = 0
    for (irank,(itrange,iurange,isrange)) in enumerate(all_ranges)
        nsites = 0
        for is in isrange, it in itrange, iu in iurange
            if iu <= it && (is+it+iu)%2 == parity
                nsites += 1
            end
        end
        min_nsites = (nsites<min_nsites) ? nsites : min_nsites
        max_nsites = (nsites>max_nsites) ? nsites : max_nsites
    end
    if max_nsites == 0
        println("max_nsites == 0 : $N,$nranks,$parity")
    end

    if min_nsites ==  typemax(Int64)
        println("min_nsites == $(typemax(Int64)) : $N,$nranks,$parity")
    end
    return (max_nsites - min_nsites)/max_nsites
end

"""Returns an estimate of the load imbalance among the ranks,
   in the range 0.0-1.0.
   This is a ~O(N^3) function that iterates on all the sites of interest.
 """
function get_imbalance(N,nranks,get_ranges_func,parity)
    all_ranges = [get_ranges_func(N,nranks,parity,irank) for irank in 0:(nranks-1)]
    get_imbalance_from_ranges(N,nranks,all_ranges,parity)
end



# TODO: use factorization logic to decide the decomposition nranks_s * nranks_tu

function _split_nranks_in_s_and_tu(nranks)
    [(nranks_s,div(nranks,nranks_s)) for nranks_s in 1:nranks if (nranks%nranks_s) == 0]
end


function get_ranges_stu(N,nranks,parity,rank)
    imbalance = 1.0

    for (nranks_s,nranks_tu) in _split_nranks_in_s_and_tu(nranks)
        rank_s = rank % nranks_s
        rank_tu = div(rank, nranks_s, RoundToZero)
        range_s = partitions(N,nranks_s)[1+rank_s]
        range_tu = _get_ranges_tu(N,nranks_tu,rank_tu)
        range_stu = (range_s, range_tu...)
    end
    # FAKE implementation
    if rank == 0
        1:N,1:N,1:N
    else
        1:0, 1:0, 1:0
    end
end



"""Returns the number of sites having iu<=it<=Ntu"""
function get_number_of_sites(Ntu)
    Int(Ntu*(Ntu+1)/2)
end

"""Returns the number of sites having iu<=it<=Ntu with a given parity."""
function _get_number_of_sites_eo(Ntu)
    nhalf = Int(floor(Ntu/2))
    total_elements = get_number_of_sites(Ntu)

    a = nhalf*(nhalf+1)
    b = total_elements - a

    even, odd = if Ntu%2 == 0 a,b else b,a end
    even, odd
end
end
