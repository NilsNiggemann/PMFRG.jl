Obsacc = 1e-14
include("UnitTests/UnitTests.jl")

using PMFRGSolve

include("../ext/PMFRGSolveMPIExt/test/MPITest/mpi-tests.jl")

##
@testset verbose = true "PMFRGSolve tests" begin
    test_IO()
    @testset verbose = true "Single Process tests" begin
        testOneLoopSolve(Obsacc)
        testTwoLoopSolve(Obsacc)
        testParquetSolve()
    end
    test_mpi_solve()
end
