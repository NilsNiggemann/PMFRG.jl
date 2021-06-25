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

    for iu in 1:N, it in 1:N, is in 1:N
        DVc[1,is,it,iu] = -DVb[1,it,is,iu]
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

function SolveFRG(Par::Params;method=DP5(),save_everystep = false,saveat = [1.,Par.Lam_min])
    State,setup = setupSystem(Par) #Package parameter and pre-allocate arrays 
    @unpack Lam_max,Lam_min,accuracy,MinimalOutput = Par

    save_func(State,Lam,integrator) = writeOutput(State,Lam,Par)

    saved_values = SavedValues(double,Observables)
    cb = SavingCallback(save_func, saved_values,save_everystep =true,tdir=-1)

    problem = ODEProblem(getDeriv!,State,(Lam_max,Lam_min),setup)
    @time sol = solve(problem,method,reltol = accuracy,abstol = accuracy, save_everystep = save_everystep,saveat = saveat,callback=cb,dt = -0.2*Par.Lam_max)
    if !MinimalOutput
        println(sol.destats)
    end
    return sol,saved_values
end

function writeOutput(State,Lam,Par)
    @unpack MinimalOutput,N,np_vec,T,usesymmetry = Par
    f_int,gamma,Va,Vb,Vc = State.x
    chi = getChi(State,Lam,Par)[:,1]

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

        Max,MaxPos = absmax(Va)
        println("Max Va",Tuple(MaxPos) ," = ", Max)
        Max,MaxPos = absmax(Vb)
        println("Max Vb",Tuple(MaxPos) ," = ", Max)
        Max,MaxPos = absmax(Vc)
        println("Max Vc",Tuple(MaxPos) ," = ", Max)
        
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
    return Observables(chi,gamma,f_int)
end