module UnitTests
    using PMFRG,Test,SpinFRGLattices
    include("DimerTest.jl")
    export test_DimerFRG,test_Gammab_Dimer,test_Gammaa_onsite,test_tu_symmetries,test_Vertex_Derivative

    include("TypeStability.jl")
    export test_Vertex_Derivative
end