function compare_arguments_post(Xexp, X)
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
                different_val_places = val .!= expval

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

using SpinFRGLattices.SquareLattice
function generate_test_params()
    System = getSquareLattice(4, [1.0, 0.1])
    Params(System, OneLoop(), T = 0.5, N = 5, accuracy = 1e-3)
end




using HDF5
using H5Zblosc
function h5deserialize(f::HDF5.File, root::String, id::Int)
    (
        X = PMFRGCore.BubbleType(
            read(f["$root/$id/X/a"]),
            read(f["$root/$id/X/b"]),
            read(f["$root/$id/X/c"]),
            read(f["$root/$id/X/Ta"]),
            read(f["$root/$id/X/Tb"]),
            read(f["$root/$id/X/Tc"]),
            read(f["$root/$id/X/Td"]),
        ),
        State = PMFRGCore.StateType(
            read(f["$root/$id/State/f_int"]),
            read(f["$root/$id/State/γ"]),
            PMFRGCore.VertexType(
                read(f["$root/$id/State/Γ/a"]),
                read(f["$root/$id/State/Γ/b"]),
                read(f["$root/$id/State/Γ/c"]),
            ),
        ),
        Deriv = PMFRGCore.StateType(
            read(f["$root/$id/Deriv/f_int"]),
            read(f["$root/$id/Deriv/γ"]),
            PMFRGCore.VertexType(
                read(f["$root/$id/Deriv/Γ/a"]),
                read(f["$root/$id/Deriv/Γ/b"]),
                read(f["$root/$id/Deriv/Γ/c"]),
            ),
        ),
        Lam = read(f["$root/$id/Lam"]),
    )
end

"""Function that can be used to store X,State,Deriv and Lam
(the non-trivially-constructible arguments of getXBubble!)
in a HDF5 file, so that h5deserialize can read them properly"""
function h5serialize!(
    f::HDF5.File,
    root::String,
    X::PMFRGCore.BubbleType,
    S::PMFRGCore.StateType,
    D::PMFRGCore.StateType,
    Lam::Float64,
    id::Int,
)
    # X
    f["$root/$id/X/a", blosc = 9] = X.a
    f["$root/$id/X/b", blosc = 9] = X.b
    f["$root/$id/X/c", blosc = 9] = X.c
    f["$root/$id/X/Ta", blosc = 9] = X.Ta
    f["$root/$id/X/Tb", blosc = 9] = X.Tb
    f["$root/$id/X/Tc", blosc = 9] = X.Tc
    f["$root/$id/X/Td", blosc = 9] = X.Td

    # S
    f["$root/$id/State/f_int", blosc = 9] = S.f_int
    f["$root/$id/State/γ", blosc = 9] = S.γ
    f["$root/$id/State/Γ/a", blosc = 9] = S.Γ.a
    f["$root/$id/State/Γ/b", blosc = 9] = S.Γ.b
    f["$root/$id/State/Γ/c", blosc = 9] = S.Γ.c

    # D
    f["$root/$id/Deriv/f_int", blosc = 9] = D.f_int
    f["$root/$id/Deriv/γ", blosc = 9] = D.γ
    f["$root/$id/Deriv/Γ/a", blosc = 9] = D.Γ.a
    f["$root/$id/Deriv/Γ/b", blosc = 9] = D.Γ.b
    f["$root/$id/Deriv/Γ/c", blosc = 9] = D.Γ.c


    # Lam
    f["$root/$id/Lam"] = Lam
end


"""This function has only been used
to convert the existing serialized data
into HDF5.
It is kept here only for reference"""
function convert_jld2_to_hdf5(data, filename)

    println("Saving to hdf5")
    f = h5open(filename, "w")
    try
        ncases = length(data["return_value"])
        f["Ncases"] = ncases
        for i = 1:ncases
            workspace, lam, _ = (data["arguments"])[i]
            workspace_post_exp, _, _ = (data["arguments_post"])[i]

            let X, State, Deriv
                (; X, State, Deriv) = workspace
                h5serialize!(f, "arguments", X, State, Deriv, lam, i)
            end


            let X, State, Deriv
                (; X, State, Deriv) = workspace_post_exp
                h5serialize!(f, "arguments_post", X, State, Deriv, lam, i)
            end


        end
    finally
        close(f)
    end

end
