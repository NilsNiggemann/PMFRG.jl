"""
This module provides partitioning logic
for the MPI partitioning of the computation
of the arrays a,b,c (and Ta,Tb,Tc,Td)
in PMFRG.getXBubblePartition!

The ranges in the s,t,u indices are computed.

The S - T domain is partitioned in a number of rectangles.
Along the S direction, the partitions are as equal as possible,
while along the T direction the partitions are adjusted
for best load balancing,
assuming that only the elements with iu<=it need to be computed
(also taking into account that only elements
for which is+it+iu has a given parity
must be computed).
"""
module MPI_Detail
include("./partition.jl")
include("./best_partition_triangle.jl")
using .BestPartitionTriangle: get_all_ranges_stu

export get_all_ranges_stu
end
