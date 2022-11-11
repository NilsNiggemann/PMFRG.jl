#_____Requires generating_data_TSweep.jl to be run first!_____
# To test the accuracy of the results, we can use a ward identity connecting the self energy and the four point vertex, both of which can be used to compute the local spin correlator
using PMFRGEvaluation
using CairoMakie #for plotting. You can use whatever plotting package you like of course

TSweepFile = "temp/PMFRG_SquareLattice_NLen=5_N=10_l1_testFile_TSweep.h5"
##
Results = PMResults(TSweepFile) 

##
let 
    fig = Figure()
    ax = Axis(fig[1,1],xlabel = L"T",ylabel = L"\delta \chi_{ii}")
    
    T = Results.T
    gamma = Results.gamma_TxN[:,1,:] #only one inequivalent site
    Chi00 = Results.Chi_TR[:,1] # local chi component
    N_extrapolate = 60 # extrapolate number of frequencies for matsubara sum
    err = wardIdentityviolation.(T,eachrow(gamma),Chi00,N_extrapolate)
    lines!(ax,T,err)

    fig
end