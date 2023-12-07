using Test
include("../../src/mpi/best_partition_triangle.jl")
using .BestPartitionTriangle: _get_ranges_tu,
    _get_number_of_sites_eo,
    _split_nranks_in_s_and_tu,
    get_ranges_stu,
    get_imbalance,
    _count_sites


@testset verbose=true "Tests for best triangle partition" begin
    @testset verbose=true "Simple manual cases" begin
        #                    N,nranks,rank
        @test _get_ranges_tu(10,     2,   0) == (1:7,1:7)
    end
    @testset "ranges cover all sites" begin
        for N in 1:10, nranks in 1:min(5,N)
            coverage = Array{Int,2}(undef,(N,N))
            coverage .= 0
            all_ranges = [ _get_ranges_tu(N,nranks,irank) for irank in 0:(nranks-1) ]
            for it in 1:N, iu in 1:it
                for (itrange,iurange) in all_ranges
                    if it in itrange && iu in iurange
                        coverage[it,iu] += 1
                    end
                end
                @test coverage[it,iu] == 1
            end
        end
    end

    @testset "ranges do not go out of 1:N range" begin
        for N in 1:10, nranks in 1:min(5,N)
            all_ranges = [ _get_ranges_tu(N,nranks,irank) for irank in 0:(nranks-1) ]
            for (itrange,iurange) in all_ranges
                @test 1 <= itrange.start <= N
                @test 1 <= itrange.stop <= N
                @test 1 <= iurange.start <= N
                @test 1 <= iurange.stop <= N
            end
        end
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

    @testset "all sites are covered, only once" begin
        for N in 2:10,nranks in 1:min(N,5), parity in 0:1
            covered = Array{Int64,3}(undef,(N,N,N))
            covered .= 0

            stu_ranges = [get_ranges_stu(N,nranks,parity,r) for r in 0:(nranks-1)]
            @testset for it in 1:N, iu in 1:it, is in 1:N
                if (it+iu+is)%2 == parity
                    for (itrange,iurange,isrange)  in stu_ranges
                        if (is in isrange &&
                            it in itrange &&
                            iu in iurange)
                            covered[is,it,iu] += 1
                        end
                    end
                    @test covered[is,it,iu] == 1
                end
            end
        end
    end

    function count_sites_reference(itrange,iurange,isrange,parity)
        nsites = 0
        for is in isrange, it in itrange, iu in iurange
            if iu <= it && (is+it+iu)%2 == parity
                nsites += 1
            end
        end
        nsites
    end

    @testset "count_sites is same as count_sites_reference" begin
        for args in [
            ( 1:1 , 1:1 ,  1:1 , 0),
            ( 1:1 , 1:1 ,  1:1 , 1),
            ( 1:3 , 1:3 ,  1:3 , 0),
            ( 1:3 , 1:3 ,  1:3 , 1),
            ( 2:4 , 1:4 ,  2:4 , 0),
            ( 2:4 , 1:4 ,  2:4 , 1),
            ( 1:35, 1:35,  1:50, 1),
            (20:35, 1:35,  1:26, 1),
            ( 1:35, 1:35,  1:25, 1),
            (26:50, 1:50,  1:25, 1),
            ( 1:35, 1:35, 26:50, 1),
            (26:50, 1:50, 26:50, 1),
            ]
            fast = _count_sites(args...)
            slow = count_sites_reference(args...)
            @test slow == fast
        end
    end

    @testset "imbalance must be zero for a singe rank" begin
        @testset for N in 2:20, parity in 0:1
            nranks = 1
            imbalance =  get_imbalance(N,nranks,get_ranges_stu,parity)
            @test imbalance == 0
        end
    end

    @testset "imbalance is <15% for relevant use cases" begin
        @testset for N in 10:5:50,nranks in 2:min(div(N,4),5), parity in 0:1
            imbalance =  get_imbalance(N,nranks,get_ranges_stu,parity)
            if imbalance >= 0.15
                println("N:$N, nranks:$nranks, parity:$parity - imbalance = $imbalance")
            end
            @test imbalance <= 0.15
        end
    end

end
