using Test
include("./decompose.jl")
using .Decompose: decompose


function test_decompose()
    @testset verbose=true "Decompose" begin
        @test decompose(6,2) == [[1,6], [2,3],[3,2],[6,1]]
        @test decompose(8,3) == [ [1,1,8],
                                  [1,2,4],
                                  [1,4,2],
                                  [1,8,1],
                                  [2,1,4],
                                  [2,2,2],
                                  [2,4,1],
                                  [4,1,2],
                                  [4,2,1],
                                  [8,1,1],
                                  ]
    end
end
