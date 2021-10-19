Base.show(io::IO, f::Float64) = @printf(io, "%1.15f", f)
##

## ______________ Helper for symmetry check ______________

function getDeriv!(Deriv,State,XandPar,Lam)
    X,XTilde,Par = XandPar #use pre-allocated X and XTilde to reduce garbage collector time
    N = Par.N
    Workspace = Workspace_Struct(Deriv,State,X,XTilde)
    @unpack DVb,DVc = Workspace
    getDFint!(Workspace,Lam,Par)
    get_Self_Energy!(Workspace,Lam,Par)
    getVertexDeriv!(Workspace,Lam,Par)

    for iu in 1:N, it in 1:N, is in 1:N, R in Par.OnsitePairs
        DVc[R,is,it,iu] = -DVb[R,it,is,iu]
    end
    return
end

function CreateX(Number,VDims)
    return ArrayPartition(Tuple(zeros(double,VDims) for _ in 1:Number))
end

function setupSystem(Par::Params)
    @unpack N,Ngamma,Npairs,VDims,couplings,T,NUnique = Par
    println("T= ",T)

    ##Allocate Memory:
    State = ArrayPartition(
        zeros(double,NUnique), # f_int 
        zeros(double,NUnique,Ngamma), # gamma
        zeros(double,VDims), #Va
        zeros(double,VDims), #Vb
        zeros(double,VDims) #Vc
    )

    X = CreateX(3,VDims)
    XTilde = CreateX(4,VDims)

    Vc = State.x[5]
    for is in 1:N, it in 1:N, iu in 1:N, Rj in 1:Npairs
        Vc[Rj,is,it,iu] = -couplings[Rj]
    end
    return State,(X,XTilde,Par)
end

function SolveFRG(Par::Params,MainFile = nothing, CheckpointDirectory = nothing;Group =string(Par.T)::String, kwargs...)
    typeof(CheckpointDirectory)==String && (CheckpointDirectory = setupDirectory(CheckpointDirectory,Par))
    State,setup = setupSystem(Par) #Package parameter and pre-allocate arrays 
    

    sol,saved_values = launchPMFRG!(State,setup,getDeriv!,Par, CheckpointDirectory = CheckpointDirectory; kwargs...)
    if MainFile !== nothing
        saveMainOutput(MainFile,sol,saved_values,Par,Group)
    end
    return sol,saved_values
end

function launchPMFRG!(State,setup,Deriv!::Function,Par::Params;CheckpointDirectory = nothing,method = DP5(),MaxVal = 50,ObsSaveat = nothing,VertexCheckpoints = [],kwargs...)
    @unpack Lam_max,Lam_min,accuracy,MinimalOutput = Par
    save_func(State,Lam,integrator) = getObservables(State,Lam,Par)
    saved_values = SavedValues(double,Observables)

    function output_func(State,Lam,integrator)
        setCheckpoint(CheckpointDirectory,State,saved_values,Lam,Par,VertexCheckpoints)
        println("") 
        writeOutput(State,saved_values,Lam,Par)
    end

    sort!(VertexCheckpoints)
    #get Default for lambda range for observables
    ObsSaveat = getLambdaMesh(ObsSaveat,Lam_min,Lam_max)

    saveCB = SavingCallback(save_func, saved_values,save_everystep =false,saveat = ObsSaveat,tdir=-1)
    outputCB = FunctionCallingCallback(output_func,tdir=-1,func_start = false)
    unstable_check(dt,u,p,t) = maximum(abs,u) >MaxVal # returns true -> Interrupts ODE integration if vertex gets too big

    problem = ODEProblem(Deriv!,State,(Lam_max,Lam_min),setup)
    #Solve ODE. default arguments may be added to, or overwritten by specifying kwargs
    @time sol = solve(problem,method,reltol = accuracy,abstol = accuracy, save_everystep = false,callback=CallbackSet(saveCB,outputCB),dt=0.2*Lam_max,dtmin = 0.1*Lam_min,unstable_check = unstable_check;kwargs...)
    if !MinimalOutput
        println(sol.destats)
    end
    return sol,saved_values
end

function getObservables(State,Lam,Par)
    @unpack MinimalOutput,N,np_vec,T,usesymmetry = Par
    f_int,gamma,Va,Vb,Vc = State.x
    chi = getChi(State,Lam,Par)
    MaxVa = @view maximum(abs,Va,dims = (2,3,4,5))[:,1,1,1]
    MaxVb = @view maximum(abs,Vb,dims = (2,3,4,5))[:,1,1,1]
    MaxVc = @view maximum(abs,Vc,dims = (2,3,4,5))[:,1,1,1]
    return Observables(chi,copy(gamma),copy(f_int),MaxVa,MaxVb,MaxVc) # make sure to allocate new memory each time this function is called
end

function writeOutput(State,saved_values,Lam,Par)
    @unpack MinimalOutput,N,np_vec,T,usesymmetry = Par
    f_int,gamma,Va,Vb,Vc = State.x
    chi = saved_values.saveval[end].Chi
    if !MinimalOutput 
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
        flush(stdout)
    end
    return 
end

function setCheckpoint(Directory::String,State,saved_values,Lam,Par,checkPointList)
    println("Time taken for output saving: ")
    @time begin
        if !isempty(checkPointList)
            if Lam < last(checkPointList) 
                Checkpoint = pop!(checkPointList)
                CheckpointFile = UniqueDirName("$Directory/$(strd(Lam)).h5")
                println("\nsaving Checkpoint Lam ≤ $Checkpoint at ",Lam)
                println("in file ", CheckpointFile)
                println("")
                mv(joinpath(Directory,"CurrentState.h5"),CheckpointFile)
            end
        end
        saveCurrentState(Directory,State,saved_values,Lam,Par)
    end
end

function setCheckpoint(Directory::Nothing,State,saved_values,Lam,Par,checkPointList)
    return
end

function getLambdaMesh(Saveat::Nothing,Lam_min,Lam_max)
    dense_range = collect(LinRange(Lam_min,5.,100))
    medium_range = collect(LinRange(5.,10.,50))
    sparse_range = collect(LinRange(10.,Lam_max,30))
    ObsSaveat = unique!(append!(dense_range,medium_range,sparse_range))
    return ObsSaveat
end

function getLambdaMesh(Saveat::AbstractVector)
    return Saveat
end