using Test

@testset verbose = true "PMFRG" begin
    include("Doctests/doctests.jl")
    @testset verbose = true "PMFRGCore" begin
        include("../PMFRGCore.jl/test/runtests.jl")
    end
    @testset verbose = true "PMFRGSolve" begin
        include("../PMFRGSolve.jl/test/runtests.jl")
    end
end
