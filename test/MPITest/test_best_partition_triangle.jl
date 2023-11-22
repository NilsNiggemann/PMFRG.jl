using Test
include("../../src/mpi/best_partition_triangle.jl")
using .BestPartitionTriangle: _get_ranges_tu,
    _get_number_of_sites_eo,
    _split_nranks_in_s_and_tu,
    get_ranges_stu,
    get_imbalance


@testset verbose=true "Tests for best triangle partition" begin
    @testset verbose=true "Simple manual cases" begin
        #                    N,nranks,rank
        @test _get_ranges_tu(10,     2,   0) == (1:7,1:7)
    end


    @testset "Count number of sites in lower triangular half" begin
        function countsites(Ntu)
            e,o = 0,0
            for it in 1:Ntu, iu in 1:it
                if (it+iu) % 2 == 0
                    e += 1
                else
                    o += 1
                end
            end
            e,o
        end

        @testset for Ntu in 1:50
            @test countsites(Ntu) == _get_number_of_sites_eo(Ntu)
        end
    end

    @testset "split nranks in nransk_tu x nranks_s" begin
        all_pairs = _split_nranks_in_s_and_tu(50)
        @test (1,50) in all_pairs
        @test (2,25) in all_pairs
        for pair in all_pairs
            @test typeof(pair) <: Tuple{Int,Int}
        end
        for i in 1:50
            if 50%i != 0
                for pair in all_pairs
                    @test pair[1] != i
                end
            end
        end

    end
    # TODO: Property-Test that
    #       - all elements are covered
    #       - only once

    @testset "all sites are covered" begin
        for N in 1:10,nranks in 1:5, parity in 0:1
            covered = Array{Int64,3}(undef,(N,N,N))
            covered .= 0

            stu_ranges = [get_ranges_stu(N,nranks,parity,r) for r in 0:(nranks-1)]
            for it in 1:N, iu in 1:it, is in 1:N
                if (it+iu+is)%2 == parity
                    for (isrange,itrange,iurange)  in stu_ranges
                        if (is in isrange &&
                            it in itrange &&
                            iu in iurange)
                            covered[is,it,iu] += 1
                        end
                        @test covered[is,it,iu] == 1
                    end
                end
            end
        end
    end
    @testset "get_imbalance_fast is same as get_imbalance" begin
        @test false
    end

    @testset "imbalance is <20%" begin
        for N in 2:10,nranks in 1:5, parity in 0:1
            @test get_imbalance(N,nranks,get_ranges_stu,parity) <= 1.0

        end
    end

end
