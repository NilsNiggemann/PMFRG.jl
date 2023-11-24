module MPI_Detail
include("./partition.jl")
include("./decompose.jl") # TODO: remove
include("./best_partition.jl") # TODO: remove
include("./best_partition_triangle.jl")
using .BestPartitionTriangle: get_all_ranges_stu

export get_all_ranges_stu
end
