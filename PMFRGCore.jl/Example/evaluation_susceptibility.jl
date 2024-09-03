#_____Requires generating_data.jl to be run first!_____
using HDF5, PMFRGEvaluation
using CairoMakie #for plotting. You can use whatever plotting package you like of course

mainFile = "temp/PMFRG_SquareLattice_NLen=5_N=10_l1_testFile.h5"
System = SquareLattice.getSquareLattice(5)
Lattice = LatticeInfo(System, SquareLattice)
let
    chi_ΛR = h5read(mainFile, "0.5/Chi")
    chi_R = chi_ΛR[:, end]

    chi = getNaiveLatticeFT(chi_R, Lattice)
    # or: chi = getLatticeFFT(chi_R, Lattice)

    k = LinRange(-2pi, 2pi, 300)

    chik = [chi(x, y) for x in k, y in k]

    fig, ax, hm = heatmap(k, k, chik, axis = (; aspect = 1))
    Colorbar(fig[1, 2], hm)
    fig
end
