using MPI
using Test

include("./test_partition.jl")
include("./test_best_partition_triangle.jl")
function test_mpi_core()
    @testset verbose = true "MPI tests" begin
        @testset verbose = true "Unit tests for MPI functionality" begin
            test_1D_partition()
            test_best_partition_triangle()
        end

        function run_mpi_script(script, n, testname)
            function print_header()
                linelength = 79
                println("="^linelength)
                title = " $testname "
                subtitle = " \"$(basename(script))\" "
                println("="^5 * title * "="^(linelength - 5 - length(title)))
                println("="^5 * subtitle * "="^(linelength - 5 - length(subtitle)))
                println("="^linelength)
            end

            print_header()
            @testset verbose = true "$testname" begin
                p = run(ignorestatus(`$(mpiexec())
                                      -n $n
                                      $(Base.julia_cmd())
                                      --project=$(Base.active_project())
                                      $script`))
                @test success(p)
            end
        end

        @testset verbose = true "MPI tests - external executables" begin
            dir = dirname(@__FILE__)
            run_mpi_script(
                joinpath(dir, "..", "RegressionTests", "PMFRG.getXBubbleMPI.jl"),
                2,
                "Regression test - getXBubbleMPI",
            )

            run_mpi_script(
                joinpath(dir, "test_chunk_communication.jl"),
                4,
                "Ibcast! communication example - test_chunk_communication.jl",
            )

        end
    end
end
