Obsacc = 1e-14
include("UnitTests/UnitTests.jl")

using MPI
@assert !isnothing(Base.get_extension(PMFRGCore, :PMFRGCoreMPIExt)) "Perhaps you need `using MPI`?"
include("../ext/PMFRGCoreMPIExt/test/MPITest/mpi-tests.jl")
include("RegressionTests/PMFRGCore.getXBubble.jl")

@testset verbose = true "PMFRGCore tests" begin
    test_mpi()
    test_getXBubble()
    testOneLoop(Obsacc)
    testTwoLoop(Obsacc)
    testParquet()
end
