Base.show(io::IO, f::Float64) = @printf(io, "%1.15f", f)
##
_getFloatType(Par::PMFRGParams) = typeof(Par.NumericalParams.T)

function InitializeState(Par::PMFRGParams)
    (; N, Ngamma) = Par.NumericalParams
    VDims = getVDims(Par)
    (; couplings, NUnique) = Par.System

    floattype = _getFloatType(Par)

    State = ArrayPartition( #Allocate Memory:
        zeros(floattype, NUnique), # f_int 
        zeros(floattype, NUnique, Ngamma), # gamma
        zeros(floattype, VDims), #Va
        zeros(floattype, VDims), #Vb
        zeros(floattype, VDims), #Vc
    )

    Γc = State.x[5]
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

    Buffs = BufferType(PropsBuffers, VertexBuffers)
    return (; X, Buffs, Par, ParallelizationScheme)
end

"""Converts t step used for integrator to Λ. Inverse of Lam_to_t."""
t_to_Lam(t) = exp(t)
"""Converts physical cutoff Λ to t (integrator step). Inverse of t_to_Lam."""
Lam_to_t(t) = log(t)

"""
Given a set of Parameters (currently, oneß and twoloop are implemented), solves the set of differential flow equations and returns the ODE solution object along with an array of Observables at each lambda step.
Allowed keyword arguments (with default values):

    MainFile = nothing,                             # Specifies name of main output file as a string.
                                                    # Defaults to 'nothing', in which case no output file is produced.
    Group = DefaultGroup(Par),               # Specifies the name of the subgroup of the main files HDF5 group to which the output data is written. Defaults to temperature
                                                    # Defaults to the value of temperature to allow temperature sweeps to be written to the same file.
    CheckpointDirectory = nothing,                  # Directory in which the current ODE state is written in regular intervals during the integration. 
                                                    # Default 'nothing' will not produce any backup!
    method = DP5(),                                 # ODE integration method. Standard is set to OrdinaryDiffEq.DP5(). 
                                                    # See OrdinaryDiffEq package for further options.
    MaxVal = Inf                                    # Terminates the ODE solution, when any absolute value of the Solution reaches MaxVal.
    ObsSaveat = nothing,                            # Specifies a vector of Λ values at which Observables are computed.
    VertexCheckpoints = [],                         # Specifies a vector of Λ values at which vertices shall be saved permanently without being overwritten. 
                                                    # Empty means only the current state will be saved. Requires CheckpointDir ≠ nothing !
    overwrite_Checkpoints = false::Bool,            # Specifies whether CheckpointDirectory is to be overwritten if it exists. 
                                                    # Defaults to false, in which case an unused name is generated.
    CheckPointSteps = 1,                            # Number of skipped steps before Checkpoint data is saved to reduce time spent on IO operations.
    kwargs...                                       # Additional keyword arguments are passed to OrdinaryDiffEq.solve.
                                                    # See the OrdinaryDiffEq documentation for further details.

"""
SolveFRG(
    Par,
    ParallelizationScheme::AbstractParallelizationScheme = MultiThreaded();
    kwargs...,
) = launchPMFRG!(
    InitializeState(Par),
    AllocateSetup(Par, ParallelizationScheme),
    getDeriv!;
    kwargs...,
)

function _get_output_callback(
    CheckPointSteps,
    CheckpointDirectory,
    saved_values,
    Par,
    VertexCheckpoints,
)
    # count number of outputs = number of steps.
    # CheckPointSteps gives the intervals in which checkpoints should be saved.
    i = 0
    function bareOutput(State, t, integrator)
        @timeit_debug "bareOutput" begin
            Lam = t_to_Lam(t)
            i += 1
            i % CheckPointSteps == 0 && setCheckpoint(
                CheckpointDirectory,
                State,
                saved_values,
                Lam,
                Par,
                VertexCheckpoints,
            )
        end
    end

    function verboseOutput(State, t, integrator)
        @timeit_debug "verboseOutput" begin
            Lam = t_to_Lam(t)
            println("Time taken for output saving: ")
            bareOutput(State, t, integrator)
            println("")
            writeOutput(State, saved_values, Lam, Par)
        end
    end

    function getOutputfunction(MinimalOutput)
        if MinimalOutput
            return bareOutput
        else
            return verboseOutput
        end
    end
    output_func = getOutputfunction(Par.Options.MinimalOutput)
    outputCB = FunctionCallingCallback(output_func, tdir = -1, func_start = false)
    return outputCB
end


function _get_saving_callback(
    ObservableType,
    Par,
    ObsSaveat,
    Lam_min,
    Lam_max,
    saved_values,
)
    #get Default for lambda range for observables
    # ObsSaveat = getLambdaMesh(ObsSaveat,Lam_min,Lam_max)
    save_func(State, t, integrator) =
        getObservables(ObservableType, State, t_to_Lam(t), Par)
    ObsSaveat = gettMesh(ObsSaveat, Lam_min, Lam_max)
    saveCB = SavingCallback(
        save_func,
        saved_values,
        save_everystep = false,
        saveat = ObsSaveat,
        tdir = -1,
    )
    return saveCB
end

function _get_problem(Lam_max, Lam_min, State, setup, Deriv!)
    t0 = Lam_to_t(Lam_max)
    tend = get_t_min(Lam_min)
    Deriv_subst! = generateSubstituteDeriv(Deriv!)
    problem = ODEProblem(Deriv_subst!, State, (t0, tend), setup)
    return problem
end


function launchPMFRG!(
    State,
    setup,
    Deriv!::Function;
    MainFile = nothing,
    Group = DefaultGroup(setup.Par),
    CheckpointDirectory::Union{String,Nothing} = nothing,
    method = DP5(),
    MaxVal = Inf,
    ObsSaveat = nothing,
    VertexCheckpoints = [],
    overwrite_Checkpoints = false::Bool,
    CheckPointSteps = 1,
    ObservableType = Observables,
    kwargs...,
)

    sort!(VertexCheckpoints)

    CheckpointDirectory =
        setupDirectory(CheckpointDirectory, setup.Par, overwrite = overwrite_Checkpoints)

    (; Lam_max, Lam_min, accuracy) = setup.Par.NumericalParams
    saved_values = SavedValues(eltype(State), ObservableType)

    #Solve ODE. default arguments may be added to, or overwritten by specifying kwargs
    setup.Par.Options.MinimalOutput || println("Starting solve")
    @timeit_debug "total solver" sol = solve(
        _get_problem(Lam_max, Lam_min, State, setup, Deriv!),
        method,
        reltol = accuracy,
        abstol = accuracy,
        save_everystep = false,
        callback = CallbackSet(
            _get_saving_callback(
                ObservableType,
                setup.Par,
                ObsSaveat,
                Lam_min,
                Lam_max,
                saved_values,
            ),
            _get_output_callback(
                CheckPointSteps,
                CheckpointDirectory,
                saved_values,
                setup.Par,
                VertexCheckpoints,
            ),
        ),
        dt = Lam_to_t(0.2 * Lam_max),
        unstable_check = (dt, u, p, t) -> maximum(abs, u) > MaxVal, # returns true -> Interrupts ODE integration if vertex gets too big
        kwargs...,
    )
    if !setup.Par.Options.MinimalOutput
        println(sol.stats)
    end
    saved_values.t .= t_to_Lam.(saved_values.t)
    saveCurrentState(
        CheckpointDirectory,
        sol.u[end],
        saved_values,
        t_to_Lam(sol.t[end]),
        setup.Par,
    )
    saveMainOutput(MainFile, sol, saved_values, setup.Par, Group)

    SetCompletionCheckmark(CheckpointDirectory)
    return sol, saved_values
end

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

function getObservables(::Type{Observables}, State::ArrayPartition, Lam, Par)
    @timeit_debug "get_observables" begin
        f_int, gamma, Va, Vb, Vc = State.x
        chi = getChi(State, Lam, Par)
        MaxVa = maximum(abs, Va, dims = (2, 3, 4, 5))[:, 1, 1, 1]
        MaxVb = maximum(abs, Vb, dims = (2, 3, 4, 5))[:, 1, 1, 1]
        MaxVc = maximum(abs, Vc, dims = (2, 3, 4, 5))[:, 1, 1, 1]
        return Observables(chi, copy(gamma), copy(f_int), MaxVa, MaxVb, MaxVc) # make sure to allocate new memory each time this function is called
    end
end

writeOutput(State::ArrayPartition, saved_values, Lam, Par) =
    writeOutput(State.x..., saved_values.saveval[end], Lam, Par)

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
