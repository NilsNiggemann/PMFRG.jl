using SpinFRGLattices, PMFRG

Par = Params( #create a group of all parameters to pass them to the FRG Solver
    getPolymer(2), # geometry, this is alw#ays required
    OneLoop(), # method. OneLoop() is the default
    N = 24, # Number of positive Matsubara frequencies for the four-point vertex.
    lenIntw = 24,
    lenIntw_acc = 24,
    accuracy = 1e-6, #absolute and relative tolerance of the ODE solver.
    # For further optional arguments, see documentation of 'NumericalParams'
    T_min = 0.5,
    T_max = exp(5),
)
struct NewObs{T}
    Chi::Vector{T}
    Chinu::Matrix{T}
    gamma::Matrix{T}
    f_int::Vector{T}
    MaxVa::Vector{T}
    MaxVb::Vector{T}
    MaxVc::Vector{T}
end

function PMFRG.getObservables(::Type{NewObs},State::PMFRG.ArrayPartition,T,Par)
    f_int,gamma,Va,Vb,Vc = State.x
    chinu = PMFRG.getChi(State,T,Par,Par.NumericalParams.N)
    MaxVa = maximum(abs,Va,dims = (2,3,4,5))[:,1,1,1]
    MaxVb = maximum(abs,Vb,dims = (2,3,4,5))[:,1,1,1]
    MaxVc = maximum(abs,Vc,dims = (2,3,4,5))[:,1,1,1]
    return NewObs(chinu[:,1],chinu,copy(gamma),copy(f_int),MaxVa,MaxVb,MaxVc) # make sure to allocate new memory each time this function is called
end

##
MainFile = nothing
# MainFile = "TFlowDimer.h5"
Solution,saved_values = SolveFRG(Par;MainFile,method = DP5(),
ObservableType = NewObs, ObsSaveat = exp10.(LinRange(log10(Par.NumericalParams.T_min),log10(100),500)))
Res = PMFRG.StructArray(saved_values.saveval) |> reverse
T = saved_values.t |> reverse
##
f_int = mean.(Res.f_int)

@. f_int = f_int*T^(0.5)

plotThermo(T,f_int,(:fint,),axkwargs = (;xscale = log10),)


##
using CairoMakie, PMFRGDimerBenchmark,PMFRGEvaluation
chi = Res.Chi
plotDimerSusc(T,chi,axkwargs = (;xscale = log10,yscale = Makie.pseudolog10))
let 
    St = PMFRG.InitializeState(Par)
    chitriv = [PMFRG.getChi(St,t, Par,1)[:,1] for t in T]
    lines!(T,getindex.(chitriv,1))
    lines!(T,getindex.(chitriv,2))
    current_figure()
end 
##
mean(x) = sum(x)/length(x)
e = let 
    Chinu = PMFRG.convertToArray(Res.Chinu)
    e = [EnergyBetaDimer(Chinu[:,:,ti], 24)*T[ti] for ti in eachindex(T)]
    
end
f_int = mean.(Res.f_int)

@. f_int = f_int

plotThermo(T,f_int,(:fint,),axkwargs = (;xscale = log10),)

# lines!(T,e,color = :darkred,linestyle = :dash)
# xlims!(1,100)
# ylims!(-0.11,0)#
current_figure()

##

gamma_nT = PMFRG.convertToArray(Res.gamma)[1,:,:] 
for i in axes(gamma_nT,1)
    @. gamma_nT[i,:] .*= T^(1/2)
end
plotGammaT(T,gamma_nT,axkwargs= (;xscale = log10,yscale = identity,
),(0,1,5))
ylims!(0,0.05)
current_figure()