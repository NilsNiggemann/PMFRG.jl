Base.show(io::IO, f::Float64) = @printf(io, "%1.15f", f)
##
_getFloatType(Par::PMFRGParams) = typeof(Par.NumericalParams.T)


function InitializeState(Par::PMFRGParams, _::MultiThreaded = MultiThreaded())
    (; couplings) = Par.System

    floattype = _getFloatType(Par)

    args = getArrayGeometry(Par)

    State = createStateVector(args; floattype = floattype)
    Γc = getVc(State, args)
    setToBareVertex!(Γc, couplings)
    return State
end

function getChannel(Buffs::AbstractVector{<:T}) where {T}
    BufferChannel = Channel{T}(length(Buffs))
    for buff in Buffs
        put!(BufferChannel, buff)
    end
    return BufferChannel
end

function AllocateSetup(
    Par::AbstractOneLoopParams,
    ParallelizationScheme::AbstractParallelizationScheme = MultiThreaded(),
)
    (; Npairs, NUnique) = Par.System
    Par.Options.MinimalOutput || println("One Loop: T= ", Par.NumericalParams.T)
    ##Allocate Memory:
    X = BubbleType(Par)
    floattype = _getFloatType(Par) #get type of float, i.e. Float64
    VertexBuffers =
        getChannel([VertexBufferType(floattype, Npairs) for _ = 1:Threads.nthreads()])
    PropsBuffers = getChannel([
        MMatrix{NUnique,NUnique,floattype,NUnique * NUnique}(undef) for
        _ = 1:Threads.nthreads()
    ])

    (; Ngamma) = Par.NumericalParams
    StateBuff = StateType(NUnique, Ngamma, getVDims(Par), floattype)
    DerivBuff = StateType(NUnique, Ngamma, getVDims(Par), floattype)

    Buffs = BufferType(PropsBuffers, VertexBuffers)
    return (; X, Buffs, Par, ParallelizationScheme, StateBuff, DerivBuff)
end

"""Converts t step used for integrator to Λ. Inverse of Lam_to_t."""
t_to_Lam(t) = exp(t)
"""Converts physical cutoff Λ to t (integrator step). Inverse of t_to_Lam."""
Lam_to_t(t) = log(t)


function generateSubstituteDeriv(getDeriv!::Function)

    function DerivSubs!(Deriv, State, par, t)
        Lam = t_to_Lam(t)
        a = getDeriv!(Deriv, State, par, Lam)
        Deriv .*= Lam
        a
    end

end


function get_t_min(Lam)
    Lam < exp(-30) && @warn "Lam_min too small! Set to exp(-30) instead."
    max(Lam_to_t(Lam), -30.0)
end

DefaultGroup(Par::PMFRGParams) = strd(Par.NumericalParams.T)


function getObservables(::Type{Observables}, State::AbstractVector, Lam, Par)
    @timeit_debug "get_observables" begin
        f_int, gamma, Va, Vb, Vc = unpackStateVector(State, Par)
        chi = getChi(State, Lam, Par)
        MaxVa = maximum(abs, Va, dims = (2, 3, 4, 5))[:, 1, 1, 1]
        MaxVb = maximum(abs, Vb, dims = (2, 3, 4, 5))[:, 1, 1, 1]
        MaxVc = maximum(abs, Vc, dims = (2, 3, 4, 5))[:, 1, 1, 1]
        return Observables(chi, copy(gamma), copy(f_int), MaxVa, MaxVb, MaxVc) # make sure to allocate new memory each time this function is called
    end
end

writeOutput(State::AbstractVector, saved_values, Lam, Par) = writeOutput(
    unpackStateVector(State, getArrayGeometry(Par))...,
    saved_values.saveval[end],
    Lam,
    Par,
)

function writeOutput(f_int, gamma, Va, Vb, Vc, obs, Lam, Par)
    (; usesymmetry) = Par.Options
    (; N, np_vec, T) = Par.NumericalParams
    chi = obs.Chi
    t = Lam_to_t(Lam)
    print(
        "T= ",
        strd(T),
        " at t step: ",
        strd(t),
        ", Λ = exp(t) = ",
        strd(Lam),
        "\tchi_1 = ",
        strd(chi[1]),
        "\tchi_2 = ",
        strd(chi[2]),
        "\t f_int = (",
    )
    for f in f_int
        print(strd(f), ",")
    end
    println(")")
    function givefreqs()
        N <= 7 && return 1, 2, 3
        f1 = 1
        f2 = max(1, div(N, 2) - 3)
        f3 = max(1, N - 5)

        n1, n2, n3 = np_vec[f1], np_vec[f2], np_vec[f3]
        while (n1 + n2 + n3) % 2 == 0 && f3 > 1
            f3 -= 1
            n3 = np_vec[f3]
        end
        return f1, f2, f3
    end
    MaxVa, MaxPosVa = absmax(Va)
    MaxVb, MaxPosVb = absmax(Vb)
    MaxVc, MaxPosVc = absmax(Vc)
    println("Max Va", Tuple(MaxPosVa), " = ", MaxVa)
    println("Max Vb", Tuple(MaxPosVb), " = ", MaxVb)
    println("Max Vc", Tuple(MaxPosVc), " = ", MaxVc)

    f1, f2, f3 = givefreqs()
    println("\t_____Symmetry tests_____")
    println("\t+Va_1($f1,$f2,$f3) = ", +Va[1, f1, f2, f3])
    println("\t-Va_1($f3,$f2,$f1) = ", -Va[1, f3, f2, f1])
    println("\t+Va_1($f2,$f3,$f1) = ", +Va[1, f2, f3, f1])

    if (!usesymmetry)
        println("\t-Va_1($f1,$f3,$f2) = ", -Va[1, f1, f3, f2], "\n")
        println("\t+Va_2($f1,$f2,$f3) = ", +Va[2, f1, f2, f3])
        println("\t-Va_2($f1,$f3,$f2) = ", -Va[2, f1, f3, f2])
        println("\t+Vb_1($f1,$f2,$f3) = ", +Vb[1, f1, f2, f3])
        println("\t-Vb_1($f1,$f3,$f2) = ", -Vb[1, f1, f3, f2], "\n")

        println(
            "\t+Va_2($f1,$f2,$f3)\n\t-Vb_2($f1,$f2,$f3)\n\t+Vc_2($f1,$f3,$f2) = ",
            (+Va[2, f1, f2, f3] - Vb[2, f1, f2, f3] + Vc[2, f1, f3, f2]),
        )
        println("\t+Vc_2($f1,$f2,$f3) = ", +Vc[2, f1, f2, f3], "\n")

        println(
            "\t+Va_1($f1,$f2,$f3)\n\t-Vb_1($f1,$f2,$f3)\n\t+Vc_1($f1,$f3,$f2) = ",
            (+Va[1, f1, f2, f3] - Vb[1, f1, f2, f3] + Vc[1, f1, f3, f2]),
        )
        println("\t+Vc_1($f1,$f2,$f3) = ", +Vc[1, f1, f2, f3], "\n")
    end
end

function getLambdaMesh(Saveat::Nothing, Lam_min, Lam_max)
    dense_range = collect(LinRange(Lam_min, 5.0, 100))
    medium_range = collect(LinRange(5.0, 10.0, 50))
    sparse_range = collect(LinRange(10.0, Lam_max, 30))
    ObsSaveat = unique!(append!(dense_range, medium_range, sparse_range))
    return ObsSaveat
end

# function gettMesh(Saveat::Nothing,Lam_min,Lam_max)
#     tmin = get_t_min(Lam_min)
#     tmax = Lam_to_t(Lam_max)
#     LinRange(tmin,tmax,150)
# end
gettMesh(Saveat, Lam_min, Lam_max) = Lam_to_t.(getLambdaMesh(Saveat, Lam_min, Lam_max))

function getLambdaMesh(Saveat::Vector{Float64}, Lam_min, Lam_max)
    return unique(push!(Saveat, Lam_max)) # make sure that there is at least one element at beginning of code
end
