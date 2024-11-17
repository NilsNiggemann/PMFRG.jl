using SpinFRGLattices, PMFRG
##
"""Define new struct of observables"""
struct NewObs{T}
    Chi::Vector{T}
    Chinu::Matrix{T}
    gamma::Matrix{T}
    f_int::Vector{T}
    MaxVa::Vector{T}
    MaxVb::Vector{T}
    MaxVc::Vector{T}
end

"""dispatch on NewObs"""
function PMFRG.getObservables(::Type{NewObs}, State::AbstractVector, T, Par)
    f_int, gamma, Va, Vb, Vc = unpackStateVector(State, Par)

    chinu = PMFRG.getChi(State, T, Par, Par.NumericalParams.N)
    MaxVa = maximum(abs, Va, dims = (2, 3, 4, 5))[:, 1, 1, 1]
    MaxVb = maximum(abs, Vb, dims = (2, 3, 4, 5))[:, 1, 1, 1]
    MaxVc = maximum(abs, Vc, dims = (2, 3, 4, 5))[:, 1, 1, 1]
    return NewObs(chinu[:, 1], chinu, copy(gamma), copy(f_int), MaxVa, MaxVb, MaxVc) # make sure to allocate new memory each time this function is called
end

System = getPolymer(2)

Trange = 0.3:0.05:1.5
for T in Trange
    Par = Params(System, OneLoop(), T = T, N = 20, accuracy = 1e-5, MinimalOutput = true)

    mainFile = "temp/" * PMFRG.generateFileName(Par, "_testFile_TSweep")

    Solution, saved_values =
        SolveFRG(Par, ObservableType = NewObs, MainFile = mainFile, method = DP5())
end
