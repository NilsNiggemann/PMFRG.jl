using PMFRG
using SpinFRGLattices.SquareLattice
using Test
using Serialization
thisdir = dirname(@__FILE__)
data = deserialize(joinpath(thisdir, "PMFRG.getXBubble.data"))

include("PMFRG.getXBubble.common.jl")

""" Produces A OneLoopParams object with no meaning, just to use for dispatch"""
function fake_oneloop_pars()::OneLoopParams
    Params(getSquareLattice(1, [0.1]), OneLoop(), T = 0.5, N = 1, accuracy = 1.0)
end


function test_getXBubble()
    @testset verbose = true "Tests for getXBubble!" begin
        @testset for i = 1:length(data["return_value"])
            return_value = (data["return_value"])[i]
            arguments = (data["arguments"])[i]
            arguments_post = (data["arguments_post"])[i]
            PMFRG.getXBubble!(arguments..., fake_oneloop_pars())
            @test compare_arguments_post(arguments, arguments_post)
        end
    end
end
