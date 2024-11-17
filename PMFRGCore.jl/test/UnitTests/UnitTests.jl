using PMFRGCore, SpinFRGLattices, Test

include("StateLib.jl")
include("DimerTest.jl")

include("FunctionCompatTests.jl")

include("TypeStability.jl")

include("ParquetTest.jl")

function testOneLoopCore(Obsacc = 1e-14)
    @testset "OneLoop" verbose = true begin
        @testset "Allocations" verbose = true begin
            test_OneLoopAllocations(Params(getPolymer(2)))
            test_OneLoopAllocations(Params(SquareKagome.getSquareKagome(4)))
        end
    end
end

function testTwoLoopCore(Obsacc = 1e-14)
    @testset "TwoLoop" verbose = true begin
        @testset "Allocations" verbose = true begin
            test_TwoLoopAllocations(Params(getPolymer(2), TwoLoop()))
            test_TwoLoopAllocations(Params(SquareKagome.getSquareKagome(4), TwoLoop()))
        end
    end
end

function testParquetCore()
    @testset "Parquet" verbose = true begin
        test_SDE()
    end
end
