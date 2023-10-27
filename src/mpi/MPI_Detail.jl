module MPI_Detail
include("./partition.jl")
include("./decompose.jl")
include("./best_partition.jl")
using .BestPartition: get_ranges 

export get_ranges
end
