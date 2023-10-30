Obsacc = 1e-14
include("UnitTests/UnitTests.jl")
include("MPITest/mpi-tests.jl")
include("RegressionTests/PMFRG.getXBubble.jl")
##
@testset verbose=true "PMFRG tests" begin
    testOneLoop(Obsacc)
    testTwoLoop(Obsacc)
    testParquet()
    test_IO()
    test_mpi()
    test_getXBubble()
end
