# This script can be run via
# ~/.julia/bin/mpiexecjl --project=./Playground
#                        -n 2
#                        julia <this file>
#

using Test
using JLD2
using MPI
using PMFRG
using SpinFRGLattices.SquareLattice

thisdir = dirname(@__FILE__)
data = load_object(joinpath(thisdir, "PMFRG.getXBubble.data"))

include("../../../../test/RegressionTests/PMFRG.getXBubble.common.jl")

MPI.Init()

@testset verbose = true "Tests for getXBubble!" begin
    @testset for i = 1:length(data["return_value"])
        workspace, lam, _ = (data["arguments"])[i]
        workspace_post_exp, lam_post_exp, _ = (data["arguments_post"])[i]
        PMFRG.getXBubble!(workspace, lam, UseMPI())
        @test compare_arguments_post((workspace_post_exp, lam_post_exp), (workspace, lam))
    end
end

MPI.Finalize()
