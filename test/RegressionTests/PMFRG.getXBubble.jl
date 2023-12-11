using PMFRG
using SpinFRGLattices.SquareLattice
using Test
using JLD2
thisdir = dirname(@__FILE__)
data = load_object(joinpath(thisdir, "PMFRG.getXBubble.data"))

include("PMFRG.getXBubble.common.jl")

""" Produces A OneLoopParams object with no meaning, just to use for dispatch"""
function fake_oneloop_pars()::OneLoopParams
    Params(getSquareLattice(1, [0.1]), OneLoop(), T = 0.5, N = 1, accuracy = 1.0)
end


function test_getXBubble()
    @testset verbose = true "Tests for getXBubble!" begin
        @testset for i = 1:length(data["return_value"])
            workspace, lam, _ = (data["arguments"])[i]
            workspace_post_exp, lam_post_exp, _ = (data["arguments_post"])[i]
            PMFRG.getXBubble!(workspace, lam, fake_oneloop_pars())
            @test compare_arguments_post(
                (workspace_post_exp, lam_post_exp),
                (workspace, lam),
            )
        end
    end
end
