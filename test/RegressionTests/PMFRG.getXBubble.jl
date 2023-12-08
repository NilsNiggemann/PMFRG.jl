using PMFRG
using Test
using Serialization
thisdir = dirname(@__FILE__)
data = deserialize(joinpath(thisdir,"PMFRG.getXBubble.data"))


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
        # This might be too strict
        @assert typeof(expval) == Array{Float64,4}
        if val != expval
            absdiff = abs.(val.-expval)
            print("Test failed: $field differs: ")
            print("Absolute Difference: ", sum(absdiff), ", ")
            print("Max Difference: ", maximum(absdiff), "\n")
            different_val_places = val .!= expval
            println("First computed values: ")
            println(first(val[different_val_places],5))
            println("First expected values: ")
            println(first(expval[different_val_places],5))
            result = false
        end
    end
    result
end

function test_getXBubble()
    @testset verbose = true "Tests for getXBubble!" begin
            @testset for i = 1:length(data["return_value"])
                    return_value = (data["return_value"])[i]
                    arguments = (data["arguments"])[i]
                    arguments_post = (data["arguments_post"])[i]
                    PMFRG.getXBubble!(arguments...)
                    @test compare_arguments_post(arguments, arguments_post)
                    end
        end
end
