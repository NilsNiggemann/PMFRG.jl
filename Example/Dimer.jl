using SpinFRGLattices, PMFRG

Par = Params( #create a group of all parameters to pass them to the FRG Solver
    getPolymer(2), # geometry, this is always required
    OneLoop(), # method. OneLoop() is the default
    N = 32, # Number of positive Matsubara frequencies for the four-point vertex.
    lenIntw = 100,
    lenIntw_acc = 200,
    accuracy = 1e-6, #absolute and relative tolerance of the ODE solver.
    # For further optional arguments, see documentation of 'NumericalParams'
    T_min = 0.2,
    T_max = 100.,
)

Solution,saved_values = SolveFRG(Par,ObsSaveat = exp10.(LinRange(log10(Par.NumericalParams.T_min),log10(Par.NumericalParams.T_max),300)))

chi1_ex(B) = (exp(B) - 1 + B )/(2*(exp(B)+3))
chi2_ex(B) = -(exp(B) - 1 - B )/(2*(exp(B)+3))
chi_ex(T) = [chi1_ex(1/T),chi2_ex(1/T)]
##
using CairoMakie
let 
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], xlabel = "T", ylabel = L"\chi",
    # xscale = log10,yscale = Makie.pseudolog10
    )
    Res = PMFRG.StructArray(saved_values.saveval)
    chi = Res.Chi
    T = saved_values.t
    println(extrema(T))
    chi1 = getindex.(chi,1)
    chi2 = getindex.(chi,2)

    lines!(ax,T, chi1, color = :red)
    lines!(ax,T, chi2, color = :blue)
    lines!(ax,T, chi1_ex.(1 ./T), color = :black)
    lines!(ax,T, chi2_ex.(1 ./T), color = :black)
    xlims!(ax,0.01,3)
    fig
end
##
let 
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], xlabel = "T", ylabel = L"\chi_{err}",
    xscale = log10,yscale = log10
    )
    Res = PMFRG.StructArray(saved_values.saveval)
    chi = Res.Chi
    T = saved_values.t
    chi1 = abs.(getindex.(chi,1) .- chi1_ex.(1 ./T))
    chi2 = abs.(getindex.(chi,2) .- chi1_ex.(1 ./T))

    lines!(ax,T, chi1, color = :red)
    lines!(ax,T, chi2, color = :blue)
    lines!(ax,T, T .^- 1, color = :black, linestyle = :dash)
    xlims!(ax,1,100)
    fig
end

##
using PMFRGEvaluation

fex(T) = -T/2*log(sum(exp(-En/T) for En in (-3/4,1/4,1/4,1/4)))

let 
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], xlabel = "T", ylabel = L"f_int",
    xscale = log10,yscale = Makie.pseudolog10
    )
    Res = PMFRG.StructArray(saved_values.saveval)
    mean(x) = sum(x)/length(x)
    f_int = mean.(Res.f_int)
    T = saved_values.t
    perm = sortperm(T)
    T = T[perm]
    f_int = f_int[perm]

    f_int = f_int*T.^(1/4)
    Thermo = PMFRGEvaluation.getThermoIntPol(T,f_int  ,1)

    fint_ex =  fex.(T) .+T*log(2)

    ThermoEx = PMFRGEvaluation.getThermoIntPol(T,fint_ex,1)
    # lines!(ax,T, f_int  , color = :black)
    lines!(ax,T, Thermo.c.(T), color = :red)
    lines!(ax,T, ThermoEx.c.(T), color = :black)
    xlims!(ax,0.01,100)
    fig
end
