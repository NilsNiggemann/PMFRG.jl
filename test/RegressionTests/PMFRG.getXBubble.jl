using PMFRG
using SpinFRGLattices.SquareLattice
using Test
using JLD2
thisdir = dirname(@__FILE__)

include("PMFRG.getXBubble.common.jl")

function test_getXBubble()
    data = load_object(joinpath(thisdir, "PMFRG.getXBubble.data"))
    @testset verbose = true "Tests for getXBubble!" begin
        @testset for i = 1:length(data["return_value"])
            workspace, lam, _ = (data["arguments"])[i]

            (;X,State,Deriv) = workspace
            Par = generate_test_params()

            (;Buffs) = PMFRG.AllocateSetup(Par)
            workspace_post_exp, _, _ = (data["arguments_post"])[i]
            PMFRG.getXBubble!(X,
                              State,
                              Deriv,
                              Par,
                              Buffs,
                              lam,
                              PMFRG.MultiThreaded())
            @test compare_arguments_post(workspace_post_exp.X, X)
        end
    end
end
