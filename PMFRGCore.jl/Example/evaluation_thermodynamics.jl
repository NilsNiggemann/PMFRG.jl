#_____Requires generating_data_TSweep.jl to be run first!_____

using PMFRGEvaluation
using CairoMakie #for plotting. You can use whatever plotting package you like of course

TSweepFile = "temp/PMFRG_SquareLattice_NLen=5_N=10_l1_testFile_TSweep.h5"

Results = PMResults(TSweepFile) #Creates a useful struct that immediately computes thermal observables. Works only if the file has the correct format!
let
    fig = Figure()
    ax_e =
        Axis(fig[1, 1], ylabel = L"E/N", xlabelvisible = false, xticklabelsvisible = false)
    ax_s =
        Axis(fig[2, 1], ylabel = L"S/N", xlabelvisible = false, xticklabelsvisible = false)
    ax_c = Axis(fig[3, 1], xlabel = L"T", ylabel = L"C/N")

    for (ax, obs) in zip((ax_e, ax_s, ax_c), (Results.e, Results.s, Results.c))
        scatterlines!(ax, Results.T, obs)
    end

    rowgap!(fig.layout, 0.1)

    fig
end
