using PMFRG
using SpinFRGLattices.SquareLattice
using Test
using HDF5
thisdir = dirname(@__FILE__)

include("PMFRG.getXBubble.common.jl")


function test_getXBubble()
    fname = joinpath(thisdir, "PMFRG.getXBubble.data.h5")
    h5file = h5open(fname, "r")
    try
        @testset verbose = true "Tests for getXBubble!" begin
            ncases = read(h5file["Ncases"])
            @testset for i = 1:ncases

                (; X, State, Deriv, Lam) = h5deserialize(h5file, "arguments", i)
                X0 = X


                sumX = sum(sum(abs.(getfield(X, field))) for field in fieldnames(typeof(X)))
                @assert sumX == 0.0 "X start value should be null"

                Par = generate_test_params()

                (; Buffs) = PMFRG.AllocateSetup(Par)
                Workspace = PMFRG.OneLoopWorkspace(State, Deriv, X0, Buffs, Par)
                PMFRG.getXBubble!(Workspace, Lam, PMFRG.MultiThreaded())

                (; X) = h5deserialize(h5file, "arguments_post", i)
                @test compare_arguments_post(X, X0)
            end
        end
    finally
        close(h5file)
    end
end
