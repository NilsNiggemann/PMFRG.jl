Obsacc = 1e-14
include("UnitTests/UnitTests.jl")

using MPI, PencilArrays
@assert !isnothing(Base.get_extension(PMFRGCore, :PMFRGCoreMPIExt)) "Perhaps you need `using MPI; using PencilArrays`?"
include("../ext/PMFRGCoreMPIExt/test/MPITest/mpi-tests.jl")
include("RegressionTests/PMFRGCore.getXBubble.jl")

@testset verbose = true "PMFRGCore tests" begin
    smoketest_function_compatibilities()
    testStateUnpacking()
    test_mpi_core()
    test_getXBubble()
    testOneLoopCore(Obsacc)
    testTwoLoopCore(Obsacc)
    testParquetCore()
end
