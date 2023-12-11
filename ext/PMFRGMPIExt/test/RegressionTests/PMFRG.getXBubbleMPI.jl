# This script can be run via
# ~/.julia/bin/mpiexecjl --project=./Playground
#                        -n 2
#                        julia <this file>
#

using Test
using Serialization
using MPI
using PMFRG
using SpinFRGLattices.SquareLattice

thisdir = dirname(@__FILE__)
data = deserialize(joinpath(thisdir, "PMFRG.getXBubble.data"))

include("../../../../test/RegressionTests/PMFRG.getXBubble.common.jl")

MPI.Init()

""" Produces A MPIOneLoopParams object with no meaning, just to use for dispatch"""
function fake_mpioneloop_pars()
    Params(getSquareLattice(1, [0.1]), UseMPI(), T = 0.5, N = 1, accuracy = 1.0)
end


@testset verbose = true "Tests for getXBubble!" begin
    @testset for i = 1:length(data["return_value"])
        return_value = (data["return_value"])[i]
        arguments = (data["arguments"])[i]
        arguments_post = (data["arguments_post"])[i]
        PMFRG.getXBubble!(arguments..., fake_mpioneloop_pars())
        @test compare_arguments_post(arguments, arguments_post)
    end
end

MPI.Finalize()
