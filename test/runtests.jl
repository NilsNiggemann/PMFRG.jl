Obsacc = 1e-14
include("UnitTests/UnitTests.jl")

using MPI
@assert !isnothing(Base.get_extension(PMFRG, :PMFRGMPIExt))
include("../ext/PMFRGMPIExt/test/MPITest/mpi-tests.jl")
include("RegressionTests/PMFRG.getXBubble.jl")

##
@testset verbose = true "PMFRG tests" begin
    test_mpi()
    testOneLoop(Obsacc)
    testTwoLoop(Obsacc)
    testParquet()
    test_IO()
    test_getXBubble()
end
