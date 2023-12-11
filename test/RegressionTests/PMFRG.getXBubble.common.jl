function extractX(args)
    Workspace, Lam = args
    # We do not care about Par because it is immutable and constant
    # and about Buffer because we do not care about the contents
    (; State, Deriv, X) = Workspace

    #(;a,b,c,Ta,Tb,Tc,Td) = X

    return X
end


function compare_arguments_post(args_post_exp, args_post)
    Xexp = args_post_exp[1].X
    X = args_post[1].X

    result = true
    for field in fieldnames(typeof(X))
        val = getfield(X, field)
        expval = getfield(Xexp, field)
        # This might be too strict
        @assert typeof(expval) == Array{Float64,4}
        if val != expval
            abs_diff = abs.(val .- expval)
            total_diff = sum(abs_diff)
            if total_diff > 1e-14
                print("Test failed: $field differs significantly: ")
                print("Absolute Difference: ", total_diff, ", ")
                print("Maximum Difference: ", maximum(abs_diff), ", ")

                println("First computed values: ")
                println(first(val[different_val_places], 5))
                println("First computed values - reinterpreted: ")
                println(reinterpret(UInt64, first(val[different_val_places], 5)))

                println("First expected values: ")
                println(first(expval[different_val_places], 5))
                println("First expected values - reinterpreted: ")
                println(reinterpret(UInt64, first(expval[different_val_places], 5)))
                result = false
            end
        end
    end
    result
end
