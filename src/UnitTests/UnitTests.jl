module UnitTests
    using PMFRG,Test,SpinFRGLattices
    include("ExampleObservables.jl")
    include("DimerTest.jl")
    export test_DimerFRG,test_Gammab_Dimer,test_Gammaa_onsite,test_tu_symmetries,test_OneLoopAllocations

    include("TypeStability.jl")
    export test_OneLoopAllocations
end