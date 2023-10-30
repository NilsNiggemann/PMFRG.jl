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
data = deserialize(joinpath(thisdir,"PMFRG.getXBubble.data"))

MPI.Init()

#"Function that flattens 'args' into a Vector{Any}, so that "
function extractX(args)
    Workspace, Lam = args
    # We do not care about Par because it is immutable and constant
    # and about Buffer because we do not care about the contents
    (;State,Deriv,X) = Workspace

    #(;a,b,c,Ta,Tb,Tc,Td) = X

    return X
end


function compare_arguments_post(args_post_exp, args_post)
    Xexp = args_post_exp[1].X
    X    = args_post[1].    X

    result = true
    for field in fieldnames(typeof(X))
        val = getfield(X,field)
        expval = getfield(Xexp,field)
        if val != expval
            print("Test failed: $field differs :")
            print("Absolute Difference : ", sum( abs.(val.-expval)), "\n")
            result = false
        end
    end
    result
end

@testset verbose = true "Tests for getXBubble!" begin
        @testset for i = 1:length(data["return_value"])
                return_value = (data["return_value"])[i]
                arguments = (data["arguments"])[i]
                arguments_post = (data["arguments_post"])[i]
                PMFRG.getXBubble!(arguments...)
                @test compare_arguments_post(arguments, arguments_post)
                end
    end

MPI.Finalize()
