Obsacc = 1e-14
include("UnitTests/UnitTests.jl")

using MPI
@assert !isnothing(Base.get_extension(PMFRG, :PMFRGMPIExt)) "Perhaps you need `using MPI`?"
include("../ext/PMFRGMPIExt/test/MPITest/mpi-tests.jl")
include("RegressionTests/PMFRG.getXBubble.jl")

##
@testset verbose = true "PMFRG tests" begin
    test_mpi()
    test_getXBubble()
    testOneLoop(Obsacc)
    testTwoLoop(Obsacc)
    testParquet()
    test_IO()
end
