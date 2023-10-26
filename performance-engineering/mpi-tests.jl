# Although this is actually recommended by
# https://juliaparallel.org/MPI.jl/stable/usage/#Writing-MPI-tests,
# It actually does not work because one needs to load the project
# when launching the julia script.
# So what I am using instead at the moment is:
#
# mpiexecjl -prepend-rank --project=./Playground -n 2 julia ./PMFRG.jl/performance-engineering/generate_data_example_mpi.jl
#
# where Playground is a project/environment where PMFRG is added as a dev dependency.

using MPI
using Test

dir=dirname(@__FILE__)

@testset "hello" begin
    n = 2  # number of processes
    script = joinpath(dir,"generate_data_example_mpi.jl")
    mpiexec() do exe  # MPI wrapper
        run(`$exe -n $n $(Base.julia_cmd()) $script`)
        # alternatively:
        # p = run(ignorestatus(`...`))
        # @test success(p)
    end
end
