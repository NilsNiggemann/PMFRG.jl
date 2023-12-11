# This script can be run via
# ~/.julia/bin/mpiexecjl --project=./Playground
#                        -n 2
#                        julia <this file>
#

using PMFRG
using Test
using Serialization
using MPI
thisdir = dirname(@__FILE__)
data = deserialize(joinpath(thisdir, "PMFRG.getXBubble.data"))

include("PMFRG.getXBubble.common.jl")

MPI.Init()

test_getXBubble()

MPI.Finalize()
