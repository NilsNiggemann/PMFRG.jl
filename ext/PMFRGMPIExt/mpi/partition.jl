
"""
Small module to compute the partitions of a range 1:N
into npieces chunks.
(See tests for properties and examples).
"""
module Partition
fences(N, npieces) = [round(Int, i * N / npieces) for i = 0:npieces]
starts(N, npieces) = 1 .+ fences(N, npieces)[1:end-1]
ends(N, npieces) = fences(N, npieces)[2:end]
partitions(N, npieces) = [s:e for (s, e) in zip(starts(N, npieces), ends(N, npieces))]
end
