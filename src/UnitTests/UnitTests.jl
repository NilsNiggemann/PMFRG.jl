# module UnitTests
using Test
# using PMFRG,SpinFRGLattices,Parameters,Test
BechmarkingParams(Method) = Params(getPolymer(2),Method,N=24,Ngamma = 24,T=0.5,accuracy = 1e-3,usesymmetry = false,Lam_min = 0.,MinimalOutput = true,lenIntw = 60)

BechmarkingParams(Method::Parquet) = Params(getPolymer(2),Method,N=24,Ngamma = 24,T=0.8,accuracy = 1e-3,usesymmetry = false,Lam_min = 0.,MinimalOutput = true,lenIntw = 60)


Params(getPolymer(2),Parquet(),N=16,T=0.7,lenIntw = 120)


include("ExampleObservables.jl")
include("DimerTest.jl")
export test_DimerFRG,test_Gammab_Dimer,test_Gammaa_onsite,test_tu_symmetries,test_DimerParquet

include("BubbleTest.jl")
export test_BubbleSymmetries

include("TypeStability.jl")
export test_OneLoopAllocations
# end