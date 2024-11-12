using PMFRGCore, SpinFRGLattices, Test

BenchmarkingParams(Method, System = getPolymer(2)) = Params(
    System,
    Method,
    T = 0.5,
    N = 10,
    Ngamma = 10,
    accuracy = 1e-3,
    Lam_min = exp(-30),
    Lam_max = 100.0,
    usesymmetry = false,
    MinimalOutput = true,
    lenIntw = 60,
    lenIntw_acc = 60,
)

BenchmarkingParams(Method::Parquet, System = getPolymer(2)) = Params(
    System,
    Method,
    T = 0.8,
    N = 24,
    Ngamma = 24,
    accuracy = 1e-3,
    Lam_min = 0.0,
    Lam_max = 100.0,
    usesymmetry = false,
    MinimalOutput = true,
    lenIntw = 60,
    lenIntw_acc = 60,
)


include("ExampleObservables.jl")
include("DimerTest.jl")
include("IOTests.jl")


function testOneLoopSolve(Obsacc = 1e-14)
    @testset "OneLoop" verbose = true begin
        @testset "Dimer" verbose = true begin
            test_DimerFRG(Obsacc = Obsacc)
        end
        @testset "Squagome" verbose = true begin
            test_SquagomeFRG(OneLoop(), Obsacc = Obsacc, tol = 1e-8)
        end
    end
end

function testTwoLoopSolve(Obsacc = 1e-14)
    @testset "TwoLoop" verbose = true begin
        @testset "Dimer" verbose = true begin
            test_DimerFRG(TwoLoop(), Obsacc = Obsacc, tol = 1e-8) # accuracy of symmetries is finite, given by length of Matsubara sum
        end
        @testset "Squagome" verbose = true begin
            test_SquagomeFRG(TwoLoop(), Obsacc = Obsacc, tol = 1e-8)
        end
    end
end

function testParquetSolve()
    @testset "Parquet" verbose = true begin
        test_DimerParquet(tol = 1e-6)
    end
end
