Base.show(io::IO, f::Float64) = @printf(io, "%1.15f", f)
##
function getDeriv!(Deriv,State,setup,Lam)
    (X,Buffs,Par) = setup #use pre-allocated X and XTilde to reduce garbage collector time
    Workspace = OneLoopWorkspace(Deriv,State,X,Buffs,Par)

    getDFint!(Workspace,Lam)
    get_Self_Energy!(Workspace,Lam)
    getVertexDeriv!(Workspace,Lam)

    symmetrizeBubble!(Workspace.X,Par)
    symmetrizeVertex!(Workspace.State.Γ,Par)
    flush(stdout)
    return
end

function getDerivVerbose!(Deriv,State,setup,Lam)
    (X,Buffs,Par) = setup #use pre-allocated X and XTilde to reduce garbage collector time
    print("Workspace:\n\t") 
    @time Workspace = OneLoopWorkspace(Deriv,State,X,Buffs,Par)
    print("getDFint:\n\t") 
    @time getDFint!(Workspace,Lam)
    print("get_Self_Energy:\n\t") 
    @time get_Self_Energy!(Workspace,Lam)
    print("getVertexDeriv:\n\t") 
    @time getVertexDeriv!(Workspace,Lam)
    print("Symmetry:\n\t") 
    @time symmetrizeX!(Workspace)
    flush(stdout)
    return
end

function InitializeState(Par::OneLoopParams)
    @unpack N,Ngamma = Par.NumericalParams
    VDims = getVDims(Par)
    @unpack couplings,NUnique = Par.System
    
    State = ArrayPartition( #Allocate Memory:
        zeros(double,NUnique), # f_int 
        zeros(double,NUnique,Ngamma), # gamma
        zeros(double,VDims), #Va
        zeros(double,VDims), #Vb
        zeros(double,VDims) #Vc
    )

    Γc = State.x[5]
    setToBareVertex!(Γc,couplings)
    return State
end

function AllocateSetup(Par::OneLoopParams)
    @unpack N,Ngamma = Par.NumericalParams

    @unpack couplings,NUnique,Npairs = Par.System
    println("One Loop: T= ",Par.NumericalParams.T)
    ##Allocate Memory:
    X = BubbleType(Par)
    PropsBuffers = [Matrix{double}(undef,NUnique,NUnique) for _ in 1:Threads.nthreads()] 
	VertexBuffers = [VertexBufferType(Npairs) for _ in 1:Threads.nthreads()]
    Buffs = BufferType(PropsBuffers,VertexBuffers) 
    return (X,Buffs,Par)
end

SolveFRG(Par;kwargs...) = launchPMFRG!(InitializeState(Par),AllocateSetup(Par),getDeriv!; kwargs...)

function launchPMFRG!(State,setup,Deriv!::Function;
    MainFile= nothing,
    Group =string(setup[end].NumericalParams.T),
    CheckpointDirectory = nothing,
    method = DP5(),
    MaxVal = 50*maximum(abs,
    setup[end].System.couplings),
    ObsSaveat = nothing,
    VertexCheckpoints = [],
    overwrite_Checkpoints = false::Bool,
    CheckPointSteps = 1,
    kwargs...)
    
    Par = setup[end]
    typeof(CheckpointDirectory)==String && (CheckpointDirectory = setupDirectory(CheckpointDirectory,Par,overwrite = overwrite_Checkpoints))

    @unpack Lam_max,Lam_min,accuracy = Par.NumericalParams
    save_func(State,Lam,integrator) = getObservables(State,Lam,Par)
    
    saved_values = SavedValues(double,Observables)
    i=0 # count number of outputs = number of steps. CheckPointSteps gives the intervals in which checkpoints should be saved.

    function bareOutput(State,Lam)
        i+=1
        i%CheckPointSteps == 0 && setCheckpoint(CheckpointDirectory,State,saved_values,Lam,Par,VertexCheckpoints)
    end
    
    function verboseOutput(State,Lam)
        println("Time taken for output saving: ")
        @time bareOutput(State,Lam)
        println("") 
        writeOutput(State,saved_values,Lam,Par)
    end

    if Par.Options.MinimalOutput 
        output_func(State,Lam,integrator) = verboseOutput(State,Lam)
    else
        output_func(State,Lam,integrator) = bareOutput(State,Lam)
    end
    sort!(VertexCheckpoints)
    #get Default for lambda range for observables
    ObsSaveat = getLambdaMesh(ObsSaveat,Lam_min,Lam_max)
    saveCB = SavingCallback(save_func, saved_values,save_everystep =false,saveat = ObsSaveat,tdir=-1)
    outputCB = FunctionCallingCallback(output_func,tdir=-1,func_start = false)
    unstable_check(dt,u,p,t) = maximum(abs,u) >MaxVal # returns true -> Interrupts ODE integration if vertex gets too big

    problem = ODEProblem(Deriv!,State,(Lam_max,Lam_min),setup)
    #Solve ODE. default arguments may be added to, or overwritten by specifying kwargs
    @time sol = solve(problem,method,reltol = accuracy,abstol = accuracy, save_everystep = false,callback=CallbackSet(saveCB,outputCB),dt=0.2*Lam_max,unstable_check = unstable_check;kwargs...)
    if !Par.Options.MinimalOutput
        println(sol.destats)
    end
    saveCurrentState(CheckpointDirectory,sol[end],saved_values,sol.t[end],Par)
    if MainFile !== nothing
        saveMainOutput(MainFile,sol,saved_values,Par,Group)
    end
    SetCompletionCheckmark(CheckpointDirectory)
    return sol,saved_values
end

function getObservables(State::ArrayPartition,Lam,Par)
    f_int,gamma,Va,Vb,Vc = State.x
    chi = getChi(State,Lam,Par)
    MaxVa = maximum(abs,Va,dims = (2,3,4,5))[:,1,1,1]
    MaxVb = maximum(abs,Vb,dims = (2,3,4,5))[:,1,1,1]
    MaxVc = maximum(abs,Vc,dims = (2,3,4,5))[:,1,1,1]
    return Observables(chi,copy(gamma),copy(f_int),MaxVa,MaxVb,MaxVc) # make sure to allocate new memory each time this function is called
end
writeOutput(State::ArrayPartition,saved_values,Lam,Par) = writeOutput(State.x...,saved_values,Lam,Par)

function writeOutput(f_int,gamma,Va,Vb,Vc,saved_values,Lam,Par)
    @unpack MinimalOutput,usesymmetry = Par.Options
    @unpack N,np_vec,T = Par.NumericalParams
    chi = saved_values.saveval[end].Chi
    print("T= ",strd(T)," at Lambda step: ",strd(Lam),"\tchi_1 = ",strd(chi[1]),"\tchi_2 = ",strd(chi[2]),"\t f_int = (")
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

function getLambdaMesh(Saveat::Nothing,Lam_min,Lam_max)
    dense_range = collect(LinRange(Lam_min,5.,100))
    medium_range = collect(LinRange(5.,10.,50))
    sparse_range = collect(LinRange(10.,Lam_max,30))
    ObsSaveat = unique!(append!(dense_range,medium_range,sparse_range))
    return ObsSaveat
end

function getLambdaMesh(Saveat::Vector{Float64},Lam_min,Lam_max)
    return unique(push!(Saveat,Lam_max)) # make sure that there is at least one element at beginning of code
end