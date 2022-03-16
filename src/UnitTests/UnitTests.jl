# module UnitTests
using Test
# using PMFRG,SpinFRGLattices,Parameters,Test
BechmarkingParams(Method) = Params(
    getPolymer(2),
    Method,
    T=0.5,
    N = 24,
    Ngamma = 24,
    accuracy = 1e-3,
    Lam_min = 0.0,
    Lam_max = 100.0,
    usesymmetry = false,
    MinimalOutput = true,
    lenIntw = 60,
    lenIntw_acc = 24
)

BechmarkingParams(Method::Parquet) = Params(
    getPolymer(2),
    Method,
    T=0.8,
    N = 24,
    Ngamma = 24,
    accuracy = 1e-3,
    Lam_min = 0.0,
    Lam_max = 100.0,
    usesymmetry = false,
    MinimalOutput = true,
    lenIntw = 60,
    lenIntw_acc = 24
)


Params(getPolymer(2),Parquet(),N=16,T=0.7,lenIntw = 120)


include("ExampleObservables.jl")
include("DimerTest.jl")
export test_DimerFRG,test_Gammab_Dimer,test_Gammaa_onsite,test_tu_symmetries,test_DimerParquet

include("BubbleTest.jl")
export test_BubbleSymmetries

include("TypeStability.jl")
export test_OneLoopAllocations
# end
