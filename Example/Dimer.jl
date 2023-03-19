using SpinFRGLattices, PMFRG

Par = Params( #create a group of all parameters to pass them to the FRG Solver
    getPolymer(2), # geometry, this is always required
    OneLoop(), # method. OneLoop() is the default
    N = 24, # Number of positive Matsubara frequencies for the four-point vertex.
    accuracy = 1e-5, #absolute and relative tolerance of the ODE solver.
    # For further optional arguments, see documentation of 'NumericalParams'
    T_min = 0.3,
    T_max = 1000.
)

Solution,saved_values = SolveFRG(Par)

chi1_ex(B) = (exp(B) - 1 + B )/(2*(exp(B)+3))
chi2_ex(B) = -(exp(B) - 1 - B )/(2*(exp(B)+3))
chi_ex(T) = [chi1_ex(1/T),chi2_ex(1/T)]

using CairoMakie
let 
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], xlabel = "T", ylabel = "Z")
    chi = getfield.(saved_values.saveval,:Chi)
    tvals = filter(<=(2), saved_values.t)
    chi1 = getindex.(chi,1)[end-length(tvals)+1:end]
    chi2 = getindex.(chi,2)[end-length(tvals)+1:end]

    lines!(ax,tvals, chi1, color = :red)
    lines!(ax,tvals, chi2, color = :blue)
    lines!(ax,tvals, chi1_ex.(1 ./tvals), color = :black)
    lines!(ax,tvals, chi2_ex.(1 ./tvals), color = :black)
    fig
end