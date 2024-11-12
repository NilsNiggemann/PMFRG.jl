
"""
Given a set of Parameters (currently, oneß and twoloop are implemented), solves the set of differential flow equations and returns the ODE solution object along with an array of Observables at each lambda step.
Allowed keyword arguments (with default values):

    MainFile = nothing,                             # Specifies name of main output file as a string.
                                                    # Defaults to 'nothing', in which case no output file is produced.
    Group = PMFRGCore.DefaultGroup(Par),                # Specifies the name of the subgroup of the main files HDF5 group to which the output data is written. Defaults to temperature
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
    PMFRGCore.InitializeState(Par),
    PMFRGCore.AllocateSetup(Par, ParallelizationScheme),
    PMFRGCore.getDeriv!;
    kwargs...,
)

function get_output_callback(
    CheckPointSteps,
    CheckpointDirectory,
    saved_values,
    Par,
    VertexCheckpoints,
)::SciMLBase.DiscreteCallback
    # count number of outputs = number of steps.
    # CheckPointSteps gives the intervals in which checkpoints should be saved.
    i = 0
    function bareOutput(State, t, integrator)
        @timeit_debug "bareOutput" begin
            Lam = PMFRGCore.t_to_Lam(t)
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
            Lam = PMFRGCore.t_to_Lam(t)
            println("Time taken for output saving: ")
            bareOutput(State, t, integrator)
            println("")
            PMFRGCore.writeOutput(State, saved_values, Lam, Par)
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


function get_saving_callback(
    ObservableType,
    Par,
    ObsSaveat,
    Lam_min,
    Lam_max,
    saved_values,
)::SciMLBase.DiscreteCallback
    #get Default for lambda range for observables
    # ObsSaveat = getLambdaMesh(ObsSaveat,Lam_min,Lam_max)
    save_func(State, t, integrator) =
        PMFRGCore.getObservables(ObservableType, State, PMFRGCore.t_to_Lam(t), Par)
    ObsSaveat = PMFRGCore.gettMesh(ObsSaveat, Lam_min, Lam_max)
    saveCB = SavingCallback(
        save_func,
        saved_values,
        save_everystep = false,
        saveat = ObsSaveat,
        tdir = -1,
    )
    return saveCB
end

function get_problem(Lam_max, Lam_min, State, setup, Deriv!)::ODEProblem
    t0 = PMFRGCore.Lam_to_t(Lam_max)
    tend = PMFRGCore.get_t_min(Lam_min)
    Deriv_subst! = PMFRGCore.generateSubstituteDeriv(Deriv!)
    problem = ODEProblem(Deriv_subst!, State, (t0, tend), setup)
    return problem
end


function launchPMFRG!(
    State,
    setup,
    Deriv!::Function;
    MainFile = nothing,
    Group = PMFRGCore.DefaultGroup(setup.Par),
    CheckpointDirectory::Union{String,Nothing} = nothing,
    method = DP5(),
    MaxVal = Inf,
    ObsSaveat = nothing,
    VertexCheckpoints = [],
    overwrite_Checkpoints = false::Bool,
    CheckPointSteps = 1,
    ObservableType = PMFRGCore.Observables,
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
        get_problem(Lam_max, Lam_min, State, setup, Deriv!),
        method,
        reltol = accuracy,
        abstol = accuracy,
        save_everystep = false,
        callback = CallbackSet(
            get_saving_callback(
                ObservableType,
                setup.Par,
                ObsSaveat,
                Lam_min,
                Lam_max,
                saved_values,
            ),
            get_output_callback(
                CheckPointSteps,
                CheckpointDirectory,
                saved_values,
                setup.Par,
                VertexCheckpoints,
            ),
        ),
        dt = PMFRGCore.Lam_to_t(0.2 * Lam_max),
        unstable_check = (dt, u, p, t) -> maximum(abs, u) > MaxVal, # returns true -> Interrupts ODE integration if vertex gets too big
        kwargs...,
    )
    if !setup.Par.Options.MinimalOutput
        println(sol.stats)
    end
    saved_values.t .= PMFRGCore.t_to_Lam.(saved_values.t)
    saveCurrentState(
        CheckpointDirectory,
        sol.u[end],
        saved_values,
        PMFRGCore.t_to_Lam(sol.t[end]),
        setup.Par,
    )
    saveMainOutput(MainFile, sol, saved_values, setup.Par, Group)

    PMFRGCore.SetCompletionCheckmark(CheckpointDirectory)
    return sol, saved_values
end
