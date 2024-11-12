Obsacc = 1e-14
include("UnitTests/UnitTests.jl")

using MPI
@assert !isnothing(Base.get_extension(PMFRGCore, :PMFRGCoreMPIExt)) "Perhaps you need `using MPI`?"
include("../ext/PMFRGCoreMPIExt/test/MPITest/mpi-tests.jl")
include("RegressionTests/PMFRGCore.getXBubble.jl")

@testset verbose = true "PMFRGCore tests" begin
    test_setup_deriv_compatibility()
    testStateUnpacking()
    test_mpi_core()
    test_getXBubble()
    testOneLoopCore(Obsacc)
    testTwoLoopCore(Obsacc)
    testParquetCore()
end
