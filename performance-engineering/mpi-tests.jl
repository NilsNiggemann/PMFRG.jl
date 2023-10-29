using MPI
using Test


include("../src/mpi/test_decompose.jl")
include("../src/mpi/test_partition.jl")
include("../src/mpi/test_best_partition.jl")

@testset verbose=true "MPI tests" begin
    @testset verbose=true "Unit tests for MPI functionality" begin
        test_1D_partition()
        test_decompose()
        test_partitioning_tools()
    end

    function run_mpi_script(script,n,testname)
        function print_header()
            linelength=79
            println("="^linelength)
            title = " $testname "
            subtitle = " \"$(basename(script))\" "
            println("="^5 * title * "="^(linelength-5-length(title)))
            println("="^5 * subtitle * "="^(linelength-5-length(subtitle)))
            println("="^linelength)
        end

        print_header()
        @testset verbose=true "$testname" begin
            mpiexec() do exe  # MPI wrapper
                p = run(ignorestatus(`$exe
                                      -n $n
                                      $(Base.julia_cmd())
                                      --project=$(Base.active_project())
                                      $script`))
                @test success(p)
            end
        end
    end

    @testset verbose=true "MPI tests - external executables" begin
        dir=dirname(@__FILE__)
        run_mpi_script(joinpath(dir,"regression-tests","PMFRG.getXBubbleMPI.jl"),
                       2,
                       "Regression test - getXBubbleMPI")

        run_mpi_script(joinpath(dir,"..","src","mpi","test_chunk_communication.jl"),
                       4,
                       "Ibcast! communication example - test_chunk_communication.jl")

        run_mpi_script(joinpath(dir,"generate_data_example_mpi.jl"),
                       2,
                       "Generate Data Example")

    end
end
