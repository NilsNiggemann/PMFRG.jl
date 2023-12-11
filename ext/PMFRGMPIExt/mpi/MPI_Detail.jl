module MPI_Detail
include("./partition.jl")
include("./best_partition_triangle.jl")
using .BestPartitionTriangle: get_all_ranges_stu

export get_all_ranges_stu
end
