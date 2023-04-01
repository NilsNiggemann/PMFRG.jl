Base.show(io::IO, f::Float64) = @printf(io, "%1.15f", f)
##
_getFloatType(Par::PMFRGParams) = typeof(Par.NumericalParams.T_min)

function InitializeState(Par::PMFRGParams)
    (;N,Ngamma) = Par.NumericalParams
    VDims = getVDims(Par)
    (;couplings,NUnique) = Par.System

    floattype = _getFloatType(Par)
    
    State = ArrayPartition( #Allocate Memory:
        zeros(floattype,NUnique), # f_int 
        zeros(floattype,NUnique,Ngamma), # gamma
        zeros(floattype,VDims), #Va
        zeros(floattype,VDims), #Vb
        zeros(floattype,VDims) #Vc
    )

    Γc = State.x[5]
    setToBareVertex!(Γc,couplings)
    return State
end

function getChannel(Buffs::AbstractVector{<:T}) where T
    BufferChannel = Channel{T}(length(Buffs))
    for buff in Buffs
        put!(BufferChannel, buff)
    end
    return BufferChannel
end

function AllocateSetup(Par::OneLoopParams)
    (;Npairs,NUnique) = Par.System
    println("One Loop:")
    ##Allocate Memory:
    X = BubbleType(Par)
    floattype = _getFloatType(Par) #get type of float, i.e. Float64
	VertexBuffers = getChannel([VertexBufferType(floattype,Npairs) for _ in 1:Threads.nthreads()])
    PropsBuffers = getChannel([MMatrix{NUnique,NUnique,floattype,NUnique*NUnique}(undef) for _ in 1:Threads.nthreads()] )

    Buffs = BufferType(PropsBuffers,VertexBuffers) 
    return (X,Buffs,Par)
end

"""Converts t step used for integrator to Λ. Inverse of T_to_t."""
t_to_T(t) = exp(t)
"""Converts physical cutoff Λ to t (integrator step). Inverse of t_to_T."""
T_to_t(t) = log(t)

"""
Given a set of Parameters (currently, oneß and twoloop are implemented), solves the set of differential flow equations and returns the ODE solution object along with an array of Observables at each T step.
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
SolveFRG(Par;kwargs...) = launchPMFRG!(InitializeState(Par),AllocateSetup(Par),getDeriv!; kwargs...)

function launchPMFRG!(State,setup,Deriv!::Function;
    MainFile = nothing,
    Group = DefaultGroup(setup[end]),
    CheckpointDirectory = nothing,
    method = DP5(),
    MaxVal = Inf,
    ObsSaveat = nothing,
    VertexCheckpoints = [],
    overwrite_Checkpoints = false::Bool,
    CheckPointSteps = 1,
    ObservableType = ObservablesChi,
    kwargs...)
    
    Par = setup[end]
    typeof(CheckpointDirectory)==String && (CheckpointDirectory = setupDirectory(CheckpointDirectory,Par,overwrite = overwrite_Checkpoints))

    (;T_max,T_min,accuracy) = Par.NumericalParams
    save_func(State,t,integrator) = getObservables(ObservableType,State,t_to_T(t),Par)
    
    saved_values = SavedValues(eltype(State),ObservableType)
    i=0 # count number of outputs = number of steps. CheckPointSteps gives the intervals in which checkpoints should be saved.

    function bareOutput(State,t,integrator)
        T = t_to_T(t)
        i+=1
        i%CheckPointSteps == 0 && setCheckpoint(CheckpointDirectory,State,saved_values,T,Par,VertexCheckpoints)
    end
    
    function verboseOutput(State,t,integrator)
        T = t_to_T(t)
        println("Time taken for output saving: ")
        @time bareOutput(State,t,integrator)
        println("") 
        writeOutput(State,saved_values,T,Par)
    end

    function getOutputfunction(MinimalOutput)
        if MinimalOutput
            return bareOutput
        else
            return verboseOutput
        end
    end
    output_func = getOutputfunction(Par.Options.MinimalOutput)
    sort!(VertexCheckpoints)
    #get Default for lambda range for observables
    # ObsSaveat = getTempMesh(ObsSaveat,T_min,T_max)
    ObsSaveat = gettMesh(ObsSaveat,T_min,T_max)
    saveCB = SavingCallback(save_func, saved_values,save_everystep =false,saveat = ObsSaveat,tdir=-1)
    outputCB = FunctionCallingCallback(output_func,tdir=-1,func_start = false)
    unstable_check(dt,u,p,t) = maximum(abs,u) >MaxVal # returns true -> Interrupts ODE integration if vertex gets too big

    t0 = T_to_t(T_max)
    tend = get_t_min(T_min)
    Deriv_subst! = generateSubstituteDeriv(Deriv!)
    problem = ODEProblem(Deriv_subst!,State,(t0,tend),setup)
    #Solve ODE. default arguments may be added to, or overwritten by specifying kwargs
    @time sol = solve(problem,method,reltol = accuracy,abstol = accuracy, save_everystep = false,callback=CallbackSet(saveCB,outputCB),dt=T_to_t(0.2*T_max),unstable_check = unstable_check;kwargs...)
    if !Par.Options.MinimalOutput
        println(sol.destats)
    end
    saved_values.t .= t_to_T.(saved_values.t)
    saveCurrentState(CheckpointDirectory,sol[end],saved_values,t_to_T(sol.t[end]),Par)
    saveMainOutput(MainFile,sol,saved_values,Par,Group)

    SetCompletionCheckmark(CheckpointDirectory)
    return sol,saved_values
end

function generateSubstituteDeriv(getDeriv!::Function)
    
    function DerivSubs!(Deriv,State,par,t)
        T = t_to_T(t)
        a = getDeriv!(Deriv,State,par,T)
        Deriv .*= T
        a
    end

end

function get_t_min(T)
    T < exp(-30) && @warn "T_min too small! Set to exp(-30) instead."
    max(T_to_t(T),-30.)
end

DefaultGroup(Par::PMFRGParams) = ""

writeOutput(State::ArrayPartition,saved_values,T,Par) = writeOutput(State.x...,saved_values.saveval[end],T,Par)

function writeOutput(f_int,gamma,Va,Vb,Vc,obs,T,Par)
    (;usesymmetry) = Par.Options
    (;N,np_vec) = Par.NumericalParams
    chi = obs.Chi
    t = T_to_t(T)
    print("T= ",strd(T)," at t step: ",strd(t),"\tchi_1 = ",strd(chi[1]),"\tchi_2 = ",strd(chi[2]),"\t f_int = (")
    for f in f_int
        print(strd(f),",")
    end
    println(")")
    function givefreqs()
        f1 = 1 
        f2 = div(N,2)-3 
        f3 = N - 5
    
        n1,n2,n3 = np_vec[f1],np_vec[f2],np_vec[f3]
        while (n1+n2+n3)%2 == 0 && f3>0
            f3 -=1
            n3 = np_vec[f3]
        end
        return f1,f2,f3
    end
    MaxVa,MaxPosVa = absmax(Va)
    MaxVb,MaxPosVb = absmax(Vb)
    MaxVc,MaxPosVc = absmax(Vc)
    println("Max Va",Tuple(MaxPosVa) ," = ", MaxVa)
    println("Max Vb",Tuple(MaxPosVb) ," = ", MaxVb)
    println("Max Vc",Tuple(MaxPosVc) ," = ", MaxVc)
    
    f1,f2,f3 = givefreqs()
    println("\t_____Symmetry tests_____")
    println("\t+Va_1($f1,$f2,$f3) = ", +Va[1,f1,f2,f3])
    println("\t-Va_1($f3,$f2,$f1) = ", -Va[1,f3,f2,f1])
    println("\t+Va_1($f2,$f3,$f1) = ", +Va[1,f2,f3,f1])

    if(!usesymmetry)
        println("\t-Va_1($f1,$f3,$f2) = ", -Va[1,f1,f3,f2] ,"\n")
        println("\t+Va_2($f1,$f2,$f3) = ", +Va[2,f1,f2,f3] )
        println("\t-Va_2($f1,$f3,$f2) = ", -Va[2,f1,f3,f2] )
        println("\t+Vb_1($f1,$f2,$f3) = ", +Vb[1,f1,f2,f3] )
        println("\t-Vb_1($f1,$f3,$f2) = ", -Vb[1,f1,f3,f2] ,"\n")

        println("\t+Va_2($f1,$f2,$f3)\n\t-Vb_2($f1,$f2,$f3)\n\t+Vc_2($f1,$f3,$f2) = ",
            (+Va[2,f1,f2,f3] -Vb[2,f1,f2,f3] +Vc[2,f1,f3,f2]))
        println("\t+Vc_2($f1,$f2,$f3) = ", +Vc[2,f1,f2,f3] ,"\n")

        println("\t+Va_1($f1,$f2,$f3)\n\t-Vb_1($f1,$f2,$f3)\n\t+Vc_1($f1,$f3,$f2) = ",
            (+Va[1,f1,f2,f3] -Vb[1,f1,f2,f3] +Vc[1,f1,f3,f2]))
        println("\t+Vc_1($f1,$f2,$f3) = ", +Vc[1,f1,f2,f3] ,"\n")
    end
end

function getTempMesh(Saveat::Nothing,T_min,T_max)
    ObsSaveatLow = exp10.(range(log10(T_min),log10(10.),length=300))
    ObsSaveat = exp10.(range(log10(10.),log10(T_max),length=201))
    return unique!(append!(ObsSaveatLow,ObsSaveat))
end

# function gettMesh(Saveat::Nothing,T_min,T_max)
#     tmin = get_t_min(T_min)
#     tmax = T_to_t(T_max)
#     LinRange(tmin,tmax,150)
# end
gettMesh(Saveat,T_min,T_max) = T_to_t.(getTempMesh(Saveat,T_min,T_max))

function getTempMesh(Saveat::Vector{Float64},T_min,T_max)
    return unique(push!(Saveat,T_max)) # make sure that there is at least one element at beginning of code
end