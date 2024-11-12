using MPI
using Test

function test_mpi_solve()
    @testset verbose = true "MPI tests" begin
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
                joinpath(dir, "generate_data_example_mpi.jl"),
                2,
                "Generate Data Example",
            )

        end
    end
end
