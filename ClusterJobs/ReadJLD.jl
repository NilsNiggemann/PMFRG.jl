include("../main_Cluster.jl")
include("../../Systems/SimpleCubic.jl")
@everywhere include("../TwoLoopPMFRG.jl")
flush(stdout)
using .SimpleCubic
##

function getData(key,path = flowPath)
    SweepFiles = filter(t -> occursin(key,t), readdir(path))
    fint = []
    Trange = []
    gamma_T = []
    Chi_T = []
    Par = Params(System = getPolymer(2))
    for filename in SweepFiles
        jldopen(path*filename, "r")  do f # open read-only (default)
            Par = f["Par"]
            Lam = Par.Lam_min
            State = f["Solution"](Lam)
            f,gamma = State.x[1:2]
            Chi = getChi(State,Lam,Par)[:,1]
            push!(Trange,Par.T)
            push!(fint,f)
            push!(gamma_T,gamma)
            push!(Chi_T,Chi)
        end
    end
    N = Par.N
    NUnique = Par.NUnique
    NLen = Par.System.NLen
    perm = sortperm(Trange)
    return Dict("Trange" => Trange[perm],"fint_Tx" => fint[perm],"gamma_TxN" => gamma_T[perm], "Chi_TR" => Chi_T[perm], "N" => N, "NLen" => NLen, "NUnique" => NUnique)#,Chi[perm]
end


function convertTensor(NestedArray)
    dim1 = length(NestedArray)
    otherDims = size(first(NestedArray))
    type = typeof(NestedArray[begin][begin])
    Tens = zeros(type,dim1,otherDims...)
    println(size(Tens))
    for i in eachindex(NestedArray)
        for j in CartesianIndices(otherDims)
            Tens[i, j.I...] = NestedArray[i][j]
        end
    end
    return Tens
end

function saveData(Data,FileName)
    println("Saving Data to: ",FileName)
    for key in keys(Dat)
        Entry = Dat[key]
        if Entry isa Vector
            Entry = convertTensor(Entry)
        end
        h5write(FileName,key,Entry)
    end
end

##
Dat = getData("Cubic","/scratch/niggeni/FlowSolutions/PMFRG/")
saveData(Dat,"Cubic_NLen=4_N=32.h5")
##
