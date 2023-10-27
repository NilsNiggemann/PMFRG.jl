using Test
include("./best_partition.jl")
using .BestPartition: find_best_ND_partition, get_ranges

@testset verbose=true "Partitioning tools" begin
@testset verbose=true "Find best ND partition" begin
    @test find_best_ND_partition((10,10,10),8) == [2,2,2]
    @test find_best_ND_partition((10,10,10),27) == [3,3,3]
    @test find_best_ND_partition((10,10,10),5) == [1,1,5]
    @test find_best_ND_partition((5,10,15),6) == [1,2,3]
end

@testset verbose=true "Get ranges" begin
    @test get_ranges((10,10,10),8,7) == [ 6:10,6:10,6:10]
    @test get_ranges((10,10,10),27,13) == [ 4:7,4:7,4:7 ] # The middle cube

end

end
