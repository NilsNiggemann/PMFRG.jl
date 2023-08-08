using PMFRG
using PMFRG.SpinFRGLattices

System = getPolymer(2) # create a structure that contains all information about the geometry of the problem. 

Par = Params( #create a group of all parameters to pass them to the FRG Solver
    System, # geometry, this is always required
    OneLoop(), # method. OneLoop() is the default
    N = 15, # Number of positive Matsubara frequencies for the four-point vertex.
    accuracy = 1e-6, #absolute and relative tolerance of the ODE solver.
    # For further optional arguments, see documentation of 'NumericalParams'
)

mainFile = "temp/" * PMFRG.generateFileName(Par, "T-flow_hybrid_Lam0") # specify a file name for main Output

Solution, saved_values = SolveFRG(
    Par,
    MainFile = mainFile,
    method = DP5(),
)

##
using CairoMakie,PMFRG.HDF5
let
    fold = "temp/PMFRG_2Polymer_N=15_l1T-flow.h5"
    fnew = "temp/PMFRG_2Polymer_N=15_l1T-flow_hybrid_Lam0.h5"

    chiold = h5read(fold,"Chi")
    chinew = h5read(fnew,"Chi")
    Told = h5read(fold,"T")
    Tnew = h5read(fnew,"T")

    fig = Figure()
    ax = Axis(fig[1,1])

    lines!(Told[100:end],chiold[1,100:end],label = "old",color = :blue)
    lines!(Told[100:end],chiold[2,100:end],label = "old",color = :blue)
    
    # return Tnew,chinew
    lines!(Tnew[100:end],chinew[1,100:end],label = "new",color = :red)
    lines!(Tnew[100:end],chinew[2,100:end],label = "new",color = :red)

    xlims!(ax,0,2)
    # ylims!(ax,0,1e-2)
    axislegend(current_axis(),position = :rt,merge = true)
    current_figure()
end
##
