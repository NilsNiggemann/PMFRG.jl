
include("./doctest_tools.jl")
include("./doctest_tests.jl")

using .MDDocTestTools
using Markdown

function test_documentation()
    test_doctest_tools()
    @testset verbose = true "Tests for the documentation" begin
        text = String(read(open("../README.md", "r")))
        readme = Markdown.parse(text)

        @testset "Readme mentions expected value generation in Testing and reference values section" begin
            section_text = extract_section(readme, "Testing and reference values")
            @test !isnothing(section_text)
            @test contains(section_text, "PMFRG_REGEN_EXPECTED_RESULTS")

        end

        @testset "Readme mentions MPI parallelization" begin
            section_text = extract_section(readme, "MPI parallelization")
            @test !isnothing(section_text)
            @test contains(
                section_text,
                "ext/PMFRGMPIExt/test/MPITest/generate_data_example.mpi.jl",
            )
        end

    end
end
