using PMFRGSolve, PMFRGCore
using MPI, PencilArrays
using Test

MPI.Init()

Obsacc = 1e-14
include("../../../../test/UnitTests/UnitTests.jl")

@testset verbose = true "Multi Process tests" begin
    testOneLoopSolve(Obsacc, ParallelizationScheme = UseMPI())
end
MPI.Finalize()
