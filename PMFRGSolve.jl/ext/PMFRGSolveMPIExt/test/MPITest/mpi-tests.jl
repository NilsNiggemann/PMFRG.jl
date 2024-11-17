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
                function print_line(string)
                    println("="^5 * string * "="^(linelength - 5 - length(string)))
                end

                print_line(title)
                print_line(subtitle)
                print_line(" Output of first rank only (out of $n)")
                println("="^linelength)

            end
            function create_mpi_shell_wrapper()
                # From: https://www.open-mpi.org/community/lists/users/2012/02/18362.php
                mpi_shell_wrapper_text = """
#!/bin/bash
ARGS=\$@
if [[ -v OMPI_COMM_WORLD_RANK ]] # OpenMPI Case
then
   MY_MPI_RANK=\$OMPI_COMM_WORLD_RANK
elif [[ -v PMI_RANK ]]           # MPICH Case, probably
then
   MY_MPI_RANK=\$PMI_RANK
fi
if [[ \$MY_MPI_RANK == 0 ]]
then
  \$ARGS
else
  \$ARGS 1>/dev/null 2>/dev/null
fi
"""
                path = "./shell-wrapper-text.sh"
                open(path, "w") do io
                    write(io, mpi_shell_wrapper_text)
                end
                run(`chmod +x $path`)
                path
            end


            print_header()
            wrapper_path = create_mpi_shell_wrapper()
            @testset verbose = true "$testname" begin
                p = run(ignorestatus(`$(mpiexec())
                                      -n $n
                                      $wrapper_path
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
