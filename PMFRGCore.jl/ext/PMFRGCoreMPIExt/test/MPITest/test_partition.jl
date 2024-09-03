using Test
include("../../mpi/partition.jl")
using .Partition: fences, starts, ends, partitions

function test_1D_partition()
    @testset verbose = true "Test for 1D partitioning" begin
        @testset verbose = true "Tests for fence" begin
            @test fences(10, 3) == [0, 3, 7, 10]  # nonzero rest
            @test fences(9, 3) == [0, 3, 6, 9]    # zero rest
        end

        @testset verbose = true "Test for starts" begin
            @test starts(10, 3) == [1, 4, 8] # nonzero rest
            @test starts(9, 3) == [1, 4, 7] # zero rest
        end

        @testset verbose = true "Test for ends" begin
            @test ends(10, 3) == [3, 7, 10] # nonzero rest
            @test ends(9, 3) == [3, 6, 9] # zero rest
        end

        @testset verbose = true "Test for partitions" begin
            @test partitions(10, 3) == [1:3, 4:7, 8:10]
            @test length(partitions(2, 3)) == 3 # One partition will be empty
        end


        @testset verbose = true "Partitions have at most 2 sizes and differ by 1" begin
            set_partitions_lengths(N, pieces) =
                Set(length(collect(p)) for p in partitions(N, pieces))

            @test set_partitions_lengths(10, 3) == Set([3, 4])
            @test set_partitions_lengths(9, 3) == Set([3])
            @test set_partitions_lengths(17, 6) == Set([2, 3])
            @test set_partitions_lengths(30, 29) == Set([1, 2])
            @test set_partitions_lengths(10, 9) == Set([1, 2])
            @test set_partitions_lengths(10, 11) == Set([0, 1]) # One partition will be empty
        end
    end
end
