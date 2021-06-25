using Distributed
@everywhere using ClusterManagers,Parameters,JLD2,HDF5,OrdinaryDiffEq

try
    const global N_tasks = parse(Int, ARGS[1])# get number of tasks from slurm script
    const global N_worker = N_tasks
    println("Workers: ",N_worker)
    addprocs(SlurmManager(N_worker))
catch
end
function adaptFileName(FullPath)
    newpath = FullPath
    vstring = "(v_"
    while isfile(newpath)
        it = findfirst(vstring,newpath)
        if it !== nothing
            currindex = parse(Int,newpath[it[end]+1:end-1])
            newpath = string(newpath[1:it.start-1],vstring,currindex+1,")")
        else
            newpath = string(newpath,vstring,1,")")
        end
    end
    return newpath
end
function DistributedTSweep(Trange,System;mainPath = "",flowPath ="",add_name ="",kwargs...)
    mainPath = string(mainPath,"PMFRG_Results/",System.Name,"/")
    if flowPath == ""
        @warn "flowPath is not set. Will write secondary output to current directory $(pwd())"
    end
    flowPath = string(flowPath,System.Name,"/")
    @unpack Name,NUnique,Npairs,NLen = System
    T_points = length(Trange)
    N = kwargs[:N]
    mainFile = adaptFileName(string(mainPath,Name,"_N=",N,add_name,".h5"))
    mkpath(mainPath)
    mkpath(flowPath)
    T_Results = pmap(T -> Trun(T,System,flowPath,add_name = add_name; kwargs...),Trange) # Distributes separet T runs over workers
    obs = T_Results[begin].saveval[end]
    Ngamma = size(obs.gamma,2)
    Chi_TR = zeros(T_points,Npairs)
    gamma_TxN = zeros(T_points,NUnique,Ngamma)
    fint_Tx = zeros(T_points,NUnique)
    for i in eachindex(Trange)
        obs = T_Results[i].saveval[end]
        Chi_TR[i,:] = obs.Chi
        gamma_TxN[i,:,:] = obs.gamma
        fint_Tx[i,:] = obs.f_int
    end
    println("Saving Data to: ",mainFile)
    h5write(mainFile,"Trange",collect(Trange))
    h5write(mainFile,"Chi_TR",Chi_TR)
    h5write(mainFile,"gamma_TxN",gamma_TxN)
    h5write(mainFile,"fint_Tx",fint_Tx)
    h5write(mainFile,"N",N)
    h5write(mainFile,"NUnique",NUnique)
    h5write(mainFile,"NLen",NLen)
    return
end

@everywhere function Trun(T,System,flowPath;add_name ="", kwargs...)
    Par = Params(System = System,T=T,MinimalOutput = false; kwargs...)
    
    Solution,saved_values = SolveFRG(Par,method = pickMethod(;kwargs))
    Filename = adaptFileName(string(flowPath, System.Name,"_N=$(Par.N)","_T=", string(round(Par.T,digits =3)),add_name,".jld2"))
    obs = saved_values.saveval[end]
    @save Filename saved_values Solution Par 
    chi = obs.Chi
    println("T=",T,"\tchi: ",chi[1]," , ",chi[2])
    return saved_values
end
@everywhere function JobArray(taskid,Trange,System;flowpath = "",add_name ="",kwargs...)
    @unpack Name,NUnique,Npairs,NLen = System
    if flowPath == ""
        @warn "flowPath is not set. Will write secondary output to current directory $(pwd())"
    end
    mainPath = string(HomeDir,"/PMFRG_Results/",Name)
    N = kwargs[:N]
    mainFile = adaptFileName(string(mainPath,Name,"_N=",N,add_name,".h5"))
    mkpath(mainPath)
    mkpath(flowPath)
    T = Trange(taskid)
    T_Result =Trun(T,System,flowpath,add_name = add_name; kwargs...) 
    obs = T_Result.saveval[end]

    Chi_R = obs.Chi
    gamma_xN= obs.gamma
    fint_x= obs.f_int

    println("Saving Data to: ",mainFile)
    h5write(mainFile,"T",T)
    h5write(mainFile,"Chi_R",Chi_R)
    h5write(mainFile,"gamma_xN",gamma_xN)
    h5write(mainFile,"fint_x",fint_x)
    h5write(mainFile,"N",N)
    h5write(mainFile,"NUnique",NUnique)
    h5write(mainFile,"NLen",NLen)
    return
end
##
@everywhere function pickMethod(;kwargs...)
    method = DP5()
    if :accuracy in keys(kwargs)
        if kwargs[:accuracy] <1E-10
            method = Vern7()
        elseif kwargs[:accuracy] <1E-4
            method = DP5()
        else kwargs[:accuracy] <1E-2
            method = BS3()
        end
    end
    return method
end
##
