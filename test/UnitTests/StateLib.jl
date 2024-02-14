using PMFRG
using Random
using Test

include("../../src/StateLib.jl")
using .StateLib

function test_state_unpacking()
    @testset "state creation and unpacking" verbose = true begin
        @testset verbose = true for array_geometry in [
            (Ngamma=3, NUnique=7, VDims=(7, 11, 11, 11)),
            (Ngamma=1, NUnique=1, VDims=(1, 1, 1, 1)),
            (Ngamma=2, NUnique=4, VDims=(4, 10, 10, 10)),]

            @testset "Size of state" begin
                expected_size = (
                    array_geometry.NUnique + # f_int
                    array_geometry.NUnique * array_geometry.Ngamma + # gamma
                    prod(array_geometry.VDims) +  #Va
                    prod(array_geometry.VDims) +  #Vb
                    prod(array_geometry.VDims)    #Vc
                )


                @test size(CreateState(array_geometry;
                    floattype=Float64)) == (expected_size,)
            end


            @testset "State can be initialized correctly" begin
                begin
                    couplings = rand(array_geometry.NUnique)
                    ss = CreateState(array_geometry;
                        floattype=Float64)

                    PMFRG.setToBareVertex!(get_Vc(ss, array_geometry), couplings)

                    @test all(-x in [0.0, couplings...] for x in get_Vc(ss, array_geometry))
                    @test any(-x in couplings for x in get_Vc(ss, array_geometry))

                end
            end
            @testset "no overlap between ranges" begin
                functions = [
                    StateLib._get_f_int_range,
                    StateLib._get_gamma_range,
                    StateLib._get_Va_range,
                    StateLib._get_Vb_range,
                    StateLib._get_Vc_range,
                ]
                @testset for f1 in functions
                    @testset for f2 in functions
                        r1 = f1(array_geometry)
                        r2 = f2(array_geometry)
                        @test (f1 == f2 ||
                               r1.start <= r1.stop < r2.start <= r2.stop ||
                               r2.start <= r2.stop < r1.start <= r1.stop)
                    end
                end
            end

            @testset "get_f_int" begin
                state = CreateState(array_geometry; floattype=Float64)

                f_int = get_f_int(state, array_geometry)
                randn!(state)
                @test all(f_int .== state[StateLib._get_f_int_range(array_geometry)])
            end

            @testset "get_gamma" begin
                state = CreateState(array_geometry; floattype=Float64)


                gamma = get_gamma(state, array_geometry)
                randn!(state)

                @test size(gamma) == (array_geometry.NUnique, array_geometry.Ngamma)
                @test all(gamma .== reshape(state[StateLib._get_gamma_range(array_geometry)],
                    (array_geometry.NUnique, array_geometry.Ngamma)))

            end

            @testset "get_Va" begin
                state = CreateState(array_geometry; floattype=Float64)
                Va = get_Va(state, array_geometry)
                randn!(state)

                @test size(Va) == array_geometry.VDims
                @test all(Va .== reshape(state[StateLib._get_Va_range(array_geometry)], array_geometry.VDims))
            end



            @testset "get_Vb" begin
                state = CreateState(array_geometry; floattype=Float64)
                Vb = get_Vb(state, array_geometry)
                randn!(state)

                @test size(Vb) == array_geometry.VDims
                @test all(Vb .== reshape(state[StateLib._get_Vb_range(array_geometry)], array_geometry.VDims))
            end

            @testset "get_Vc" begin
                couplings = rand(array_geometry.NUnique)
                state = CreateState(array_geometry; floattype=Float64)

                Vc = get_Vc(state, array_geometry)
                randn!(state)

                @test size(Vc) == array_geometry.VDims
                @test all(Vc .== reshape(state[StateLib._get_Vc_range(array_geometry)], array_geometry.VDims))
            end


            @testset "get_all" verbose = true begin
                couplings = rand(array_geometry.NUnique)
                state = CreateState(array_geometry; floattype=Float64)

                all_parts = unpack_state_vector(state, array_geometry)

                fs = [get_f_int, get_gamma, get_Va, get_Vb, get_Vc]
                @testset "equivalence to single functions" for (f, part) in zip(fs, all_parts)
                    @test all(part .== f(state, array_geometry))
                end

                @testset "view, no copy" begin
                    state[1] += 354.8

                    @test all_parts[1][1] == state[1]
                end


            end
        end

    end
end
