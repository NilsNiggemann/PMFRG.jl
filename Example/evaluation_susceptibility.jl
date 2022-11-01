using HDF5, FRGLatticeEvaluation
using CairoMakie #for plotting. You can use whatever plotting package you like of course
##
Lattice = LatticeInfo(System,SquareLattice)

chi_ΛR = h5read(mainFile,"0.5/Chi")
chi_R = chi_ΛR[:,end]

chi = getFourier(chi_R,Lattice)

k = LinRange(-2pi,2pi,100)

chik = [chi(x,y) for x in k, y in k]

heatmap(k,k,chik)