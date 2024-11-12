# This script contains just an example that can be run with
# ~/.julia/bin/mpiexecjl --project=./Playground
#                        -n 4
#                        julia <this file>

using Test
using MPI

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
nranks = MPI.Comm_size(comm)

@assert nranks == 4


Npairs = 12
N = 10
full_array = fill(rank, Npairs, N, N, N)

nprange = 1:Npairs
isrange = 1:N

function ranges_per_rank(root)
    itidx = floor(Int, root / 2)
    itrange = Int(N / 2 * itidx + 1):Int(N / 2 * (itidx + 1))
    iuidx = root % 2
    iurange = Int(N / 2 * iuidx + 1):Int(N / 2 * (iuidx + 1))
    itrange, iurange
end

if rank == 0
    for other_rank = 0:3
        print("Ranges for rank $other_rank:", ranges_per_rank(other_rank), "\n")
    end
end

full_array[nprange, isrange, ranges_per_rank(rank)...] .= rank

for root = 0:3
    itrange, iurange = ranges_per_rank(root)
    MPI.Bcast!((@view full_array[nprange, isrange, itrange, iurange]), root, comm)
end

for printer = 0:3
    MPI.Barrier(comm)
    if rank == printer
        @testset verbose = true "Test from rank $printer" begin
            @testset verbose = true for root = 0:3
                itrange, iurange = ranges_per_rank(root)
                @test all(full_array[nprange, isrange, itrange, iurange] .== root)
            end
        end
    end
end


MPI.Finalize()
