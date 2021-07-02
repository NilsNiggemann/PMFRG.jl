module Evaluation

export deriv,get_e,get_c,get_e,get_c,getNumberFromName,ReadPMResults,Energy,GetThermo,reverseTOrder,GetBetaThermo,PMResults,Thermoplots,plotgamma_T,plotgamma,getHTSE,cutData,cutDataAndRecompute

using HDF5,Plots,Parameters,LaTeXStrings,SmoothingSplines,DelimitedFiles,Dierckx
function deriv(y::AbstractArray,x,order=1) 
    func = Spline1D(x, y, k=3,bc="extrapolate") 
    # func = Spline1D(x, y; w=ones(length(x)), k=3, bc="nearest", s=0.0)
    derivative(func,x,nu=order)
end

function get_e(f,T)
    return(f(T)-T*derivative(f,T))
end

function get_c(f,T)
    return(-T*derivative(f,T,nu = 2))
end

function get_e(f::AbstractArray,T::AbstractArray)
    return( f .-T .*deriv(f,T))
end

function get_c(f::AbstractArray,T::AbstractArray)
    return( -T .* deriv(f,T,2))
end

@with_kw struct PMResults
    T::Vector{Float64}
    Chi_TR::Array{Float64,2}
    gamma_TxN::Array{Float64,3}
    fint::Vector{Float64}
    N::Int
    NLen::Int
    NUnique::Int
    f::Vector{Float64}
    e::Vector{Float64}
    c::Vector{Float64}
    s::Vector{Float64}
end

function getNumberFromName(Name,subName)
    res_string = split(Name,subName*"=")[end]
    for i in length(res_string):-1:1
        N = tryparse(Int,res_string[1:i])
        if N !== nothing
            return N
        end
    end
    error("Could not get ", subName, "from string ",Name)
end

function ReadPMResults(Filename)

    T = h5read(Filename,"Trange")
    # sort values to ascending order in T
    keylist = sortperm(T)
    T = T[keylist]

    fint_T = mean(h5read(Filename,"fint_Tx")[keylist,:],dims=2)[:,1]
    Chi_TR = h5read(Filename,"Chi_TR")[keylist,:]
    gamma_TxN = h5read(Filename,"gamma_TxN")[keylist,:,:]
    N = h5read(Filename,"N")
    NLen= 0
    try
        NLen = h5read(Filename,"NLen")
    catch
        NLen = getNumberFromName(Filename,"NLen")
    end
    NUnique = h5read(Filename,"NUnique")
    # skip values 
   
    return (Dict(:T => T ,:fint_T => fint_T ,:Chi_TR => Chi_TR ,:gamma_TxN => gamma_TxN ,:N => N ,:NLen => NLen ,:NUnique => NUnique))
end

"""Compute energy from spin correlations"""
function Energy(Chi_R, Lattice)
    @unpack PairList,SiteList,PairTypes,Basis,UnitCell,pairToInequiv = Lattice
    J_ij = Lattice.System.couplings
    E = 0.
    for i_site in UnitCell
        # println(Ri)
        for j_site in SiteList # site summation
            R_Ref,ij = pairToInequiv(i_site,j_site) #Map j to correct pair so that we may use Chi_0,j'
            xi = getSiteType(R_Ref,Basis)
            pair = MapToPair(xi,ij,PairList,PairTypes)
            if pair !== 0
                E += 3/(2*Basis.NCell) *J_ij[pair] * Chi_R[pair]
            end
            # println(j_site,Chi_R[pair])
        end
    end
    return E
end

function GetThermo(PMData::Dict;skipvals = 1,smoothen = false,smoothParam = 0.001)

    T,fint_T,Chi_TR,gamma_TxN,N,NLen,NUnique = getindex.(Ref(PMData),(:T,:fint_T,:Chi_TR,:gamma_TxN,:N,:NLen,:NUnique))

    T = T[1:skipvals:end]
    fint_T = fint_T[1:skipvals:end]
    Chi_TR = Chi_TR[1:skipvals:end,:]
    gamma_TxN = gamma_TxN[1:skipvals:end,:,:]

    f_T = similar(fint_T)
    e_T = similar(fint_T)
    c_T = similar(fint_T)
    s_T = similar(fint_T)
    # smooth Data for f:
    if smoothen
        spl = fit(SmoothingSpline, T, fint_T, smoothParam) # smoothing parameter low means less smoothing
        fint_T = SmoothingSplines.predict(spl,T) # fitted vector
    end
    f = -T*log(2) +fint_T
    # f_intPol = intpol(fint_T,T)
    f_intPol = Spline1D(T, f, k=3,bc="extrapolate") 
    # return f_intPol
    for (iT,Temp) in enumerate(T)
        e_T[iT] = get_e(f_intPol,Temp)
        c_T[iT] = get_c(f_intPol,Temp)
        s_T[iT] = (e_T[iT]-f[iT])/Temp
    end
    return PMResults(T=T,N = N, NLen = NLen, NUnique = NUnique, Chi_TR=Chi_TR,gamma_TxN=gamma_TxN,fint=fint_T,f=f,e=e_T,c=c_T,s=s_T)
end

function reverseTOrder(T,Arrays...)
    neworder = length(T):-1:1
    reorder(x) = x[neworder,fill(:,ndims(x)-1)...]
    return reorder.((T,Arrays...))
end

function PMResults(Filename;kwargs...)
    res = ReadPMResults(Filename)
    return GetThermo(res;kwargs...)

end

function Thermoplots(Results,pl =plot(layout = (4,1));method = plot!,shape = :circle,kwargs...)
    @unpack T,f,e,s,c = Results
    ThermQuantities = (f,e,s,c)
    # linestyles = (:solid,:dash,:dot,:dashdot)
    # shapes = (:circle,:rect,:diamond,:cross)
    # colors = ("blue","red","black","cyan","magenta","green","pink")
    Labels = [L"f",L"e",L"s",L"c"]
    legendLabel = (kwargs[:label],"","","")
    for (i,(obs,lab)) in enumerate(zip(ThermQuantities,Labels))
        method(pl[i],T,obs,ylabel = lab,xlabel = "",shape=shape,xformatter=_->"",top_margin = -20*Plots.px,xlims =  [0,maximum(T)];kwargs...,label = legendLabel[i])
    end
    # plot!(pl[1],,legend = true;kwargs...)
    plot!(pl[end],xlabel = L"T",xformatter=x->x,size = (500,700),link=:x,left_margin=20*Plots.px;kwargs...)
    return pl
end
function plotgamma_T(Results,iT,pl = plot())
    @unpack N,gamma_TxN = Results
    for x in 1:NUnique
        scatter!(pl,1:N,gamma_TxN[iT,x,:],ylabel = L"\gamma",xlabel = L"N",label = nothing)
    end
end

function plotgamma(Results,x=1,Nmax=size(Results.gamma_TxN,3))
    @unpack T,gamma_TxN = Results
    surface(1:Nmax,T,gamma_TxN[:,x,1:Nmax],zlabel = L"\gamma",ylabel=L"T",xlabel = L"N",label = nothing,c= colorscheme)
end


HTSE_keys = (
    "T",
    "Chi_10",
    "U_9",
    "C/k_10",
    "Chi_Pade_4,6",
    "Chi_Pade_5,5",
    "Chi_Pade_6,4",
    "U_Pade_4,5",
    "U_Pade_5,4",
    "C/k_Pade_4,6",
    "C/k_Pade_5,5",
    "C/k_Pade_6,4",
    )
function getHTSE(FileName)
    HTSE_Data = readdlm(FileName,skipstart = 155)
    HTSE = Dict( key => HTSE_Data[:,i] for (i,key) in enumerate(HTSE_keys))
end

"""Removes T points and re-computes the derivative"""
function cutData(Results,index1,index2 = 0)
    @unpack T,fint,f,e,s,c,Chi_TR,N,NLen,NUnique,gamma_TxN = Results
    fields = deepcopy((T,fint,f,e,s,c,Chi_TR,N,NLen,NUnique,gamma_TxN))
    Names = (:T,:fint,:f,:e,:s,:c,:Chi_TR,:N,:NLen,:NUnique,:gamma_TxN)
    init = Dict(Name => field for (Name,field) in zip(Names,fields))

    slice(x) = x[begin+index1:end-index2,fill(:,ndims(x)-1)...] #slices array along first dim (Temperature)
    for (key,val) in init
        if val isa AbstractArray
            init[key] = slice(val)
        end
    end
    return PMResults(;init...)
end

"""Removes T points and re-computes the derivative"""
function cutDataAndRecompute(Results,index1,index2 = 0;kwargs...)
    @unpack T,fint,Chi_TR,N,NLen,NUnique,gamma_TxN = Results
    fields = Dict(:T =>T,:fint_Tx => fint,:Chi_TR => Chi_TR,:N => N,:NLen => NLen,:NUnique => NUnique,:gamma_TxN => gamma_TxN)

    slice(x) = x[begin+index1:end-index2,fill(:,ndims(x)-1)...] #slices array along first dim (Temperature)
    for (key,val) in fields
        if val isa AbstractArray
            fields[key] = slice(val)
        end
    end
    return GetThermo(fields;kwargs...)
end

end #module