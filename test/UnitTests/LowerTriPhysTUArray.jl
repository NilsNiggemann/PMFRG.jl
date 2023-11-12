using Test
include("../../src/LowerTriPhysTUArray.jl")

@testset verbose=true "LowerHalfPhysTUArray tests" begin
    @testset verbose=true "size and length" begin
        arr = LowerTriPhysTUArray{Float64}(3,4,5,Even)
        @test Base.size(arr) == (3,4,5,5)
        @test Base.length(arr) == 3*4*5*(5+1)/2 / 2

        #arr = SymmPhysTUArray{Float64}(3,3,5,Even)
    end

    @testset verbose=true "show does not throw exceptions" begin
        arr = LowerTriPhysTUArray{Float64}(3,4,5,Even)
        tmpbuf = IOBuffer()
        @test show(tmpbuf,arr) isa Any
    end

    @testset verbose=false "Length calculation" begin
	    function calc_n_e_o(Ns,Ntu)
            e,o = 0,0
            for iu in 1:Ntu, it in 1:iu, is in 1:Ns
                if (iu+it+is)%2 == 0
                    e += 1
                else
                    o += 1
                end
            end
            e,o
        end

        @testset for (Npairs,Ns,Ntu,parity) in [(3,4,5,Even),
                                                (3,4,5,Odd),
                                                (5,25,25,Even),
                                                (5,25,25,Odd),
                                                (5,24,24,Even),
                                                (5,24,24,Odd),
                                                (5,25,24,Even),
                                                (5,25,24,Odd),
                                                (5,24,25,Even),
                                                (5,24,25,Odd)]
            arr = LowerTriPhysTUArray{Float64}(Npairs,Ns,Ntu,parity)
            parity_idx = if parity == Even 1 else 2 end
            @test Base.length(arr) == calc_n_e_o(Ns,Ntu)[parity_idx] * Npairs
        end

    end

    @testset verbose=true "indexing functions" begin
        @testset verbose=true "[s,t,u] sector" begin
            @testset "number of even and odd sites" begin
                function count_sites(Ns,Ntu)
                    e,o = 0,0
                    for iu in 1:Ntu, it in 1:iu, is in 1:Ns
                        if (iu+it+is) % 2 == 0
                            e += 1
                        else
                            o += 1
                        end
                    end
                    e,o
                end

                @testset "countsites" begin
                    @testset for Ntu in 0:50
                        @testset for Ns in 0:50
                            @test count_sites(Ns,Ntu) == n_eo_elements_in_tri_half_stu(Ns,Ntu)
                        end
                    end
                end
            end

            @testset verbose=true "iu < it or throws TriHalfError" begin
                    @test_throws TriHalfError _stu_eo_offset(3,  1,  1,  2)
            end
            @testset verbose=true "eo - simple cases" begin
                @testset "manual" begin
                    @test _stu_eo_offset(3,  1,  1,  1) == 0 # o
                    @test _stu_eo_offset(3,  2,  1,  1) == 0 # e
                    @test _stu_eo_offset(3,  3,  1,  1) == 1 # o
                    @test _stu_eo_offset(3,  1,  2,  1) == 1 # e
                    @test _stu_eo_offset(3,  2,  2,  1) == 2 # o
                    @test _stu_eo_offset(3,  3,  2,  1) == 2 # e
                    @test _stu_eo_offset(3,  1,  2,  2) == 3 # o
                    @test _stu_eo_offset(3,  2,  2,  2) == 3 # e
                    @test _stu_eo_offset(3,  3,  2,  2) == 4 # o
                    @test _stu_eo_offset(3,  1,  3,  1) == 5 # o
                    @test _stu_eo_offset(3,  2,  3,  1) == 4 # e
                    @test _stu_eo_offset(3,  3,  3,  1) == 6 # o
                    @test _stu_eo_offset(3,  1,  3,  2) == 5 # e
                    @test _stu_eo_offset(3,  2,  3,  2) == 7 # o
                    @test _stu_eo_offset(3,  3,  3,  2) == 6 # e
                    @test _stu_eo_offset(3,  1,  3,  3) == 8 # o
                    @test _stu_eo_offset(3,  2,  3,  3) == 7 # e
                    @test _stu_eo_offset(3,  3,  3,  3) == 9 # o
                end
                @testset "automatic" begin
                    @testset for (s_extent,tu_extent) in [(25,7),
                                                          (7,25),
                                                          (24,7),
                                                          (7,24),
                                                          (24,20),
                                                          (20,24),
                                                          (24,24),
                                                          (25,25)]
                        idxe, idxo = 0,0
                        for it in 1:tu_extent
                            for iu in 1:it
                                for is in 1:s_extent
                                    compare_with = if (is+iu+it)%2 == 0 idxe += 1 else idxo += 1 end
                                    @test _stu_eo_offset(s_extent,is,it,iu) == compare_with-1
                                end
                            end
                        end
                    end
                end
            end
        end

        @testset verbose=true "[r,s,t,u] sector" begin
            #             rij_extent, rij, s_extent, is, it, iu
            @test _rstu_eo_offset(15,   5,        3,  2,  2,  2) == 3*15+4

            #          rij_extent, rij, s_extent, is, it, iu
            @test _rstu_eo_idx(15,   5,        3,  2,  2,  2) == 3*15+5
        end
    end

    @testset verbose=true "setindex! and getindex" begin
        @testset verbose=true "proper constraints on parity" begin
            @testset "Even type" begin
                arr = LowerTriPhysTUArray{Float64}(3,4,5,Even)
                @test_throws ParityError arr[2,3,4,4]
                @test arr[2,3,4,3] isa Any # does not throw
            end
            @testset "Odd type" begin
                arr = LowerTriPhysTUArray{Float64}(3,4,5,Odd)
                @test_throws ParityError arr[2,3,4,3]
                @test arr[2,3,4,4] isa Any # does not throw
            end
        end
        @testset verbose=true "what you put is what you get back" begin
            @testset verbose=true "Even partition" begin
                @testset "same site" begin
                    arr = LowerTriPhysTUArray{Float64}(3,4,5,Even)
                    arr[2,3,4,3] = 42.0
                    @test arr[2,3,4,3] == 42.0
                end
                @testset "does conversion" begin
                    arr = LowerTriPhysTUArray{Float64}(3,4,5,Even)
                    arr[2,3,4,3] = 42 # Int
                    @test arr[2,3,4,3] == 42.0
                end
            end
        end
    end
end
