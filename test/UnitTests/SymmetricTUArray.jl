using Test
include("../../src/SymmetricTUArray.jl")

@testset verbose=true "SymmetricTUArray tests" begin
    @testset verbose=true "size and length" begin
        arr = SymmPhysTUArray{Float64}(3,4,5,Even)
        @test Base.size(arr) == (3,4,5,5)
        @test Base.length(arr) == 3*4*5*(5+1)/2 / 2

        #arr = SymmPhysTUArray{Float64}(3,3,5,Even)
    end

    @testset verbose=true "Length calculation" begin
	    function calc_n_e_o(N)
            e,o = 0,0
            for i in 1:N, j in 1:i
                if (i+j)%2 == 0
                    e += 1
                else
                    o += 1
                end
            end
            e,o
        end
    end

    @testset verbose=true "indexing functions" begin
        # offset only for the t,u dimensions (no eo )
        #               it, iu
        @test _tu_offset(1,  1) == 0
        @test _tu_offset(2,  1) == 1
        @test _tu_offset(1,  2) == 1
        @test _tu_offset(2,  2) == 2 # (*)

        @testset "[it,iu] sector - check symmetry" begin
            @testset verbose=false "check symmetry" for it in 1:4 for iu in 1:4
                @test _tu_offset(it,iu) == _tu_offset(iu,it)
            end
            end
        end

        @testset "[it,iu] sector - check surjectivity" begin
            # All the values in 1:(N*N+1)/2 are touched
            N = 25
            v = zeros(Int32,Int(N*(N+1)/2))
            for iu in 1:N
                for it in 1:N
                    v[_tu_offset(it,iu)+1] += 1
                end
            end
            for iu in 1:N
                for it in 1:N
                   @test v[_tu_offset(it,iu)+1] == if iu == it 1 else 2 end
                end
            end
        end

        @testset "[it,iu] sector - monotonic in iu and it" begin
            N = 25
            idx = -1
            for iu in 1:N
                for it in 1:iu
                   next_idx = _tu_offset(it,iu)
                   @test next_idx > idx
                   idx = next_idx
                end
            end
        end

        @testset "[it,iu] sector - number of even and odd sites" begin
            function count_sites(N)
                e,o = 0,0
                for iu in 1:N, it in 1:iu
                    if (iu+it) % 2 == 0
                        e += 1
                    else
                        o += 1
                    end
                end
                e,o
            end

            @testset for N in 0:20
                @test count_sites(N) == n_eo_elements_in_lower_half(N)
            end
        end

        @testset verbose=true "[it,iu] sector - eo - simple cases" begin
            @testset "manual" begin
                @test _tu_eo_offset(1,1) == 0 # e
                @test _tu_eo_offset(1,2) == 0 # o
                @test _tu_eo_offset(2,2) == 1 # e
                @test _tu_eo_offset(1,3) == 2 # e
                @test _tu_eo_offset(2,3) == 1 # o
                @test _tu_eo_offset(3,3) == 3 # e
                @test _tu_eo_offset(1,4) == 2 # o
                @test _tu_eo_offset(2,4) == 4 # e
                @test _tu_eo_offset(3,4) == 3 # o
                @test _tu_eo_offset(4,4) == 5 # e
                @test _tu_eo_offset(1,5) == 6 # e
                @test _tu_eo_offset(2,5) == 4 # o
                @test _tu_eo_offset(3,5) == 7 # e
                @test _tu_eo_offset(4,5) == 5 # o
                @test _tu_eo_offset(5,5) == 8 # e
            end

            @testset "automatic" begin
                idxe, idxo = 0,0
                for iu in 1:20
                    for it in 1:iu
                        compare_with = if (iu+it)%2 == 0 idxe += 1 else idxo += 1 end
                        @test _tu_eo_offset(it,iu) == compare_with-1
                    end
                end
            end
        end

        @testset verbose=true "[it,iu] sector - eo - check surjectivity" begin
            N = 8
            ne,no = n_eo_elements_in_lower_half(N)
            @testset "even" begin
                v = zeros(Int32,ne)
                for iu in 1:N
                    for it in 1:iu
                        if (iu+it)%2 == 0
                            @test v[_tu_eo_offset(it,iu)+1] == 0
                            v[_tu_eo_offset(it,iu)+1] += 1
                        end
                    end
                end
            end
            @testset "odd" begin
                v = zeros(Int32,no)
                for iu in 1:N
                    for it in 1:iu
                        if (iu+it)%2 == 1
                            @test v[_tu_eo_offset(it,iu)+1] == 0
                            v[_tu_eo_offset(it,iu)+1] += 1
                        end
                    end
                end
            end
        end


        @testset verbose=true "[it,iu] sector - eo - monotonic in iu and it" begin
            N = 8
            @testset "even" begin
                idx = -1
                for iu in 1:N
                    for it in 1:iu
                        if (it+iu) % 2 == 0
                            next_idx = _tu_eo_offset(it,iu)
                            @test next_idx > idx
                            idx = next_idx
                        end
                    end
                end
            end
            @testset "odd" begin
                idx = -1
                for iu in 1:N
                    for it in 1:iu
                        if (it+iu) % 2 == 1
                            next_idx = _tu_eo_offset(it,iu)
                            @test next_idx > idx
                            idx = next_idx
                        end
                    end
                end
            end
        end



        #          s_extent, is, it, iu
        @test _stu_offset(3,  1,  1,  1) == 0
        @test _stu_offset(3,  2,  1,  1) == 1
        @test _stu_offset(3,  3,  1,  1) == 2
        @test _stu_offset(3,  3,  2,  2) == 2 + 2*3 # (*)

        #             s_extent, is, it, iu
        @test _stu_eo_offset(3,  2,  2,  2) == 3
        @test _stu_eo_offset(3,  3,  2,  2) == 4

        @testset verbose=true "[is,it,iu] sector - eo - simple cases" begin
            s_extent = 5
            idxe, idxo = 0,0
            for iu in 1:s_extent
                for it in 1:iu
                    for is in 1:s_extent
                        compare_with = if (is+iu+it)%2 == 0 idxe += 1 else idxo += 1 end
                    @test _stu_eo_offset(s_extent,is,it,iu) == compare_with-1
                    end
                end
            end
        end

        @testset verbose=true "[it,iu] sector - eo - check surjectivity" begin
            N = 8
            ne,no = n_eo_elements_in_lower_half(N)
            @testset "even" begin
                v = zeros(Int32,ne)
                for iu in 1:N
                    for it in 1:iu
                        if (iu+it)%2 == 0
                            @test v[_tu_eo_offset(it,iu)+1] == 0
                            v[_tu_eo_offset(it,iu)+1] += 1
                        end
                    end
                end
            end
            @testset "odd" begin
                v = zeros(Int32,no)
                for iu in 1:N
                    for it in 1:iu
                        if (iu+it)%2 == 1
                            @test v[_tu_eo_offset(it,iu)+1] == 0
                            v[_tu_eo_offset(it,iu)+1] += 1
                        end
                    end
                end
            end
        end


        @testset verbose=true "[it,iu] sector - eo - monotonic in iu and it" begin
            N = 8
            @testset "even" begin
                idx = -1
                for iu in 1:N
                    for it in 1:iu
                        if (it+iu) % 2 == 0
                            next_idx = _tu_eo_offset(it,iu)
                            @test next_idx > idx
                            idx = next_idx
                        end
                    end
                end
            end
            @testset "odd" begin
                idx = -1
                for iu in 1:N
                    for it in 1:iu
                        if (it+iu) % 2 == 1
                            next_idx = _tu_eo_offset(it,iu)
                            @test next_idx > idx
                            idx = next_idx
                        end
                    end
                end
            end
        end




        #             rij_extent, rij, s_extent, is, it, iu
        @test _rstu_eo_offset(15,   5,        3,  2,  2,  2) == 3*15+4

        #          rij_extent, rij, s_extent, is, it, iu
        @test _rstu_eo_idx(15,   5,        3,  2,  2,  2) == 3*15+5
    end

    @testset verbose=true "setindex! and getindex" begin
        @testset verbose=true "proper constraints on parity" begin
            @testset "Even type" begin
                arr = SymmPhysTUArray{Float64}(3,4,5,Even)
                @test_throws BoundsError arr[2,3,4,4]
                @test arr[2,3,4,3] isa Any
            end
            @testset "Odd type" begin
                arr = SymmPhysTUArray{Float64}(3,4,5,Odd)
                @test_throws BoundsError arr[2,3,4,3]
                @test arr[2,3,4,4] isa Any
            end
        end
        @test begin
            arr = SymmPhysTUArray{Float64}(3,4,5,Even)
            arr[2,3,4,3] = 42.0
            arr[2,3,4,3] == 42.0
        end

        @test begin
            arr = SymmPhysTUArray{Float64}(3,4,5,Even)
            arr[2,3,1,4] = 42.0
            arr[2,3,4,1] == 42.0
        end
    end
end
