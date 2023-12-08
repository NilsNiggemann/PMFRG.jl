Obsacc = 1e-14
include("UnitTests/UnitTests.jl")
if Base.get_extension(PMFRG, :PMFRGMPI) != nothing
    println("PMFRGMPI extension loaded, will run MPI tests")
    include("MPITest/mpi-tests.jl")
else
    println("PMFRGMPI extension NOT loaded, will NOT run MPI tests")
end
include("RegressionTests/PMFRG.getXBubble.jl")

##
@testset verbose = true "PMFRG tests" begin
    testOneLoop(Obsacc)
    testTwoLoop(Obsacc)
    testParquet()
    test_IO()
    if Base.get_extension(PMFRG, :PMFRGMPI) != nothing
        test_mpi()
    end
    test_getXBubble()
end
