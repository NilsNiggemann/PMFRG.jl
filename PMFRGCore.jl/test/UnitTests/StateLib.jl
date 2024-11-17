using PMFRGCore
using Random
using Test

include("../../src/StateLib.jl")
using .StateLib

function testStateUnpacking()
    @testset "state creation and unpacking" verbose = true begin
        @testset verbose = true for array_geometry in [
            (Ngamma = 3, NUnique = 7, VDims = (7, 11, 11, 11)),
            (Ngamma = 1, NUnique = 1, VDims = (1, 1, 1, 1)),
            (Ngamma = 2, NUnique = 4, VDims = (4, 10, 10, 10)),
        ]

            @testset "Size of state" begin
                expected_size = (
                    array_geometry.NUnique + # f_int
                    array_geometry.NUnique * array_geometry.Ngamma + # gamma
                    prod(array_geometry.VDims) +  #Va
                    prod(array_geometry.VDims) +  #Vb
                    prod(array_geometry.VDims)    #Vc
                )


                @test size(createStateVector(array_geometry; floattype = Float64)) ==
                      (expected_size,)
            end


            @testset "State can be initialized correctly" begin
                begin
                    couplings = rand(array_geometry.NUnique)
                    ss = createStateVector(array_geometry; floattype = Float64)

                    PMFRGCore.setToBareVertex!(getVc(ss, array_geometry), couplings)

                    @test all(-x in [0.0, couplings...] for x in getVc(ss, array_geometry))
                    @test any(-x in couplings for x in getVc(ss, array_geometry))

                end
            end
            @testset "no overlap between ranges" begin
                functions = [
                    StateLib._getF_intRange,
                    StateLib._getGammaRange,
                    StateLib._getVaRange,
                    StateLib._getVbRange,
                    StateLib._getVcRange,
                ]
                @testset for f1 in functions
                    @testset for f2 in functions
                        r1 = f1(array_geometry)
                        r2 = f2(array_geometry)
                        @test (
                            f1 == f2 ||
                            r1.start <= r1.stop < r2.start <= r2.stop ||
                            r2.start <= r2.stop < r1.start <= r1.stop
                        )
                    end
                end
            end

            @testset "getF_int" begin
                state = createStateVector(array_geometry; floattype = Float64)

                f_int = getF_int(state, array_geometry)
                for i in eachindex(state)
                    state[i] = i
                end
                @test all(f_int .== state[StateLib._getF_intRange(array_geometry)])
            end

            @testset "getGamma" begin
                state = createStateVector(array_geometry; floattype = Float64)


                gamma = getGamma(state, array_geometry)
                for i in eachindex(state)
                    state[i] = i
                end
                @test size(gamma) == (array_geometry.NUnique, array_geometry.Ngamma)
                @test all(
                    gamma .== reshape(
                        state[StateLib._getGammaRange(array_geometry)],
                        (array_geometry.NUnique, array_geometry.Ngamma),
                    ),
                )

            end

            @testset "getVa" begin
                state = createStateVector(array_geometry; floattype = Float64)
                Va = getVa(state, array_geometry)
                for i in eachindex(state)
                    state[i] = i
                end
                @test size(Va) == array_geometry.VDims
                @test all(
                    Va .== reshape(
                        state[StateLib._getVaRange(array_geometry)],
                        array_geometry.VDims,
                    ),
                )
            end



            @testset "getVb" begin
                state = createStateVector(array_geometry; floattype = Float64)
                Vb = getVb(state, array_geometry)
                for i in eachindex(state)
                    state[i] = i
                end

                @test size(Vb) == array_geometry.VDims
                @test all(
                    Vb .== reshape(
                        state[StateLib._getVbRange(array_geometry)],
                        array_geometry.VDims,
                    ),
                )
            end

            @testset "getVc" begin
                couplings = rand(array_geometry.NUnique)
                state = createStateVector(array_geometry; floattype = Float64)

                Vc = getVc(state, array_geometry)
                for i in eachindex(state)
                    state[i] = i
                end
                @test size(Vc) == array_geometry.VDims
                @test all(
                    Vc .== reshape(
                        state[StateLib._getVcRange(array_geometry)],
                        array_geometry.VDims,
                    ),
                )
            end


        end

        @testset "packing-unpacking between State vector and StateType" verbose = true begin

            array_geometry = (NUnique = 3, Ngamma = 4, VDims = (5, 7, 6, 8))

            @testset "unpackStateVector" verbose = true begin
                couplings = rand(array_geometry.NUnique)
                state = createStateVector(array_geometry; floattype = Float64)

                all_parts = unpackStateVector(state, array_geometry)

                fs = [getF_int, getGamma, getVa, getVb, getVc]
                @testset "equivalence to single functions" begin
                    @testset for (f, part) in zip(fs, all_parts)
                        @test all(part .== f(state, array_geometry))
                    end
                end

                @testset "view, no copy" begin
                    state[1] += 354.8

                    @test all_parts[1][1] == state[1]
                end

                @testset "unpack into StateType" begin
                    S = PMFRGCore.StateType(
                        array_geometry.NUnique,
                        array_geometry.Ngamma,
                        array_geometry.VDims,
                    )
                    for i in eachindex(state)
                        state[i] = i
                    end
                    unpackStateVector!(S, state)

                    @test all(S.f_int .== getF_int(state, array_geometry))
                    @test all(S.γ .== getGamma(state, array_geometry))
                    @test all(S.Γ.a .== getVa(state, array_geometry))
                    @test all(S.Γ.b .== getVb(state, array_geometry))
                    @test all(S.Γ.c .== getVc(state, array_geometry))
                end
            end

            @testset "repackStateVector and variants" verbose = true begin
                S = PMFRGCore.StateType(
                    array_geometry.NUnique,
                    array_geometry.Ngamma,
                    array_geometry.VDims,
                )

                v = 0.0
                for i in eachindex(S.f_int)
                    S.f_int[i] = v
                    v += 1
                end
                for i in eachindex(S.γ)
                    S.γ[i] = v
                    v += 1
                end
                for i in eachindex(S.Γ.a)
                    S.Γ.a[i] = v
                    v += 1
                end
                for i in eachindex(S.Γ.b)
                    S.Γ.b[i] = v
                    v += 1
                end
                for i in eachindex(S.Γ.c)
                    S.Γ.c[i] = v
                    v += 1
                end

                @testset verbose = true "repackStateVector" begin
                    packed = repackStateVector(S)
                    all_parts = unpackStateVector(packed, array_geometry)

                    fs = [getF_int, getGamma, getVa, getVb, getVc]
                    @testset "equivalence to single functions" begin
                        @testset for (f, part) in zip(fs, all_parts)
                            @test all(part .== f(packed, array_geometry))
                        end
                    end
                end

                @testset "repackStateVector, with range" begin
                    @testset for nparts = 1:5
                        packed = repackStateVector(S)
                        total_length = length(packed)
                        chunk_length = total_length / nparts

                        roundtoint(x) = Int(round(x))

                        chunks = [
                            range(
                                start = roundtoint(chunk_length * (i - 1)) + 1,
                                stop = roundtoint(chunk_length * i),
                            ) for i = 1:nparts
                        ]

                        @testset verbose = true for range in chunks
                            buff = zeros(length(range))
                            repackStateVector!(buff, range, S)

                            @test buff == packed[range]
                        end
                    end
                end

            end

            @testset "get array geometry from state" begin
                S = PMFRGCore.StateType(
                    array_geometry.NUnique,
                    array_geometry.Ngamma,
                    array_geometry.VDims,
                )

                ag = getArrayGeometry(S)

                @test ag.Ngamma == array_geometry.Ngamma
                @test ag.NUnique == array_geometry.NUnique
                @test ag.VDims == array_geometry.VDims

            end
        end
    end # @testset "state creation and unpacking" verbose = true begin
end # function testStateUnpacking()
