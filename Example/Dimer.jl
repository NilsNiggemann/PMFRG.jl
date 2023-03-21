using SpinFRGLattices, PMFRG

Par = Params( #create a group of all parameters to pass them to the FRG Solver
    getPolymer(2), # geometry, this is alw#ays required
    OneLoop(), # method. OneLoop() is the default
    N = 20, # Number of positive Matsubara frequencies for the four-point vertex.
    lenIntw = 20,
    lenIntw_acc = 20,
    accuracy = 1e-6, #absolute and relative tolerance of the ODE solver.
    # For further optional arguments, see documentation of 'NumericalParams'
    T_min = 0.05,
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

Solution,saved_values = SolveFRG(Par,MainFile = "TFlowDimer.h5",method = DP5(),
ObservableType = NewObs, ObsSaveat = exp10.(LinRange(log10(Par.NumericalParams.T_min),log10(100),500)))
Res = PMFRG.StructArray(saved_values.saveval) |> reverse
T = saved_values.t |> reverse
##
using CairoMakie, PMFRGDimerBenchmark,PMFRGEvaluation
chi = Res.Chi
plotDimerSusc(T,chi,axkwargs = (;xscale = log10,yscale = Makie.pseudolog10))
##
mean(x) = sum(x)/length(x)
e = let 
    Chinu = PMFRG.convertToArray(Res.Chinu)
    e = [get_e_Chi(Chinu)
    
end
f_int = mean.(Res.f_int)

f_int = f_int .*T.^(3/3)

plotThermo(T,f_int,axkwargs = (;xscale = log10))
##
gamma_nT = PMFRG.convertToArray(Res.gamma)[1,:,:] 
for i in axes(gamma_nT,1)
    @. gamma_nT[i,:] .*= T^(3/2)
end
plotGammaT(T,gamma_nT,axkwargs= (;xscale = log10,yscale = identity,
),(0,1,5))
current_figure()