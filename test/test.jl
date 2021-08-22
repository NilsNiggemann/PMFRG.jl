using PMFRG,Profile,SpinFRGLattices
# using SpinFRGLattices.Octochlore
using SpinFRGLattices.SimpleCubic
##
println("Start program")
flush(stdout)
##
System = getCubic(4)
function Testrun(System)
    # System = getPyrochlore(4)
    # System.couplings ./= maximum(System.couplings)
    # System = getCubic(5)
    # Par = Params(System = System,N=32,T=2.2,MinimalOutput=false,usesymmetry = true,accuracy= 1E-4,Ngamma = 200,Lam_min=145.)
    # SolveFRG(Par,method = BS3())
    # println("Setup-run includes compile time")

    Par = Params(System = System,N=24,T=1.2,MinimalOutput=false,usesymmetry = true,accuracy= 1E-4,Ngamma = 32)
    # Solution,saved_values = SolveFRG(Par,method = BS3(),CheckpointDirectory = "/scratch/niggeni/test",VertexCheckpoints = [10.5,0.2])
    Solution,saved_values = SolveFRG_Checkpoint("F:/scratch/niggeni/test/Cubic_NLen=4_N=24_T=1.200000000000000/Checkpoints/10.401.h5",System,method = BS3(),CheckpointDirectory = "/scratch/niggeni/test",VertexCheckpoints = [10.5,0.2])
    return Par,Solution,saved_values
end

##
# @profview Par,Solution,saved_values = Testrun()
Par,Solution,saved_values = Testrun(System)
# State = Solution(0.01)
# println(@Free_Energy())
obs = saved_values.saveval[end]
# println(obs.f_int)
println(obs.Chi)
##
using Plots,PMFRGEvaluation,RecursiveArrayTools,StaticArrays
Oct = getOctochlore(6,0,0)
Lattice = LatticeInfo(System = Oct,Basis = Octochlore.Basis,pairToInequiv = Octochlore.pairToInequiv)
chi = VectorOfArray(getfield.(saved_values.saveval,:Chi))
k,Jk = Fourier2D(-Oct.couplings,(x,y)-> SA[x,x,y],Lattice,ext= 4pi,res= 200)
Chikplot(k,Jk,dpi =500)
##
plotMaxFlow(chi',saved_values.t,Lattice,(x,y) -> SA[x,x,y])
##
w = PMFRG.get_w.(1:Par.Ngamma,Ref(0.2))
# scatter(saved_values.t,getindex.(getfield.(saved_values.saveval,:Chi),1))
scatter(w,saved_values.saveval[end].gamma[1,:])
plot!(w,9 ./(16 .*w),label = "\$T=0\$")
vline!([w[Par.N]],linestyle = :dash)
##
# Par = Params(System=getDimerSquareKagome(3,[1.,1,1,0]),N=64);
Par = Params(System=getCubic(10),N=20);
using BenchmarkTools,Parameters,RecursiveArrayTools
##
function test(Par)
    @unpack VDims,NUnique,Ngamma,T = Par
    State = ArrayPartition(
        zeros(double,NUnique), # f_int 
        zeros(double,NUnique,Ngamma), # gamma
        zeros(double,VDims), #Va
        zeros(double,VDims), #Vb
        zeros(double,VDims) #Vc
    )

    Deriv = similar(State);
    X =  ArrayPartition(
        zeros(double,VDims),
        zeros(double,VDims),
        zeros(double,VDims)
    )
    XTilde =  ArrayPartition(
        zeros(double,VDims),
        zeros(double,VDims),
        zeros(double,VDims),
        zeros(double,VDims)
    )
    @allocated Workspace = PMFRG.Workspace_Struct(Deriv,State,X,XTilde)
    Lam = 10.

    # display(@benchmark getDeriv!($Deriv,$State,$(X,XTilde,Par),$Lam) evals = 6)
    # @btime PMFRG.getVertexDeriv!($Workspace,$Lam,$Par)
    Buf = PMFRG.VertexBuffer(Par.Npairs)


    # @btime PMFRG.addX!($Workspace,1,1,1,2,$[1],$Par,$Buf)
    # display(@benchmark PMFRG.getVertexDeriv!($Workspace,$Lam,$Par) evals = 1)
    @btime PMFRG.getVertexDeriv!($Workspace,$Lam,$Par)
    # display(@benchmark get_Self_Energy!($Workspace,$Lam,$Par) evals = 6)
    # @benchmark getChi($State, $Lam,$Par)
    # for _ in 1:10
    #     println(@time getDeriv!(Deriv,State,(X,XTilde,Par),Lam))
    # end
    # @code_warntype sChannel!(Workspace,3,3,3,1.,0.1,Par)
    # @btime get_Self_Energy!($Workspace,0.1,$Par)

    # @benchmark getKataninProp($gamma,$Dgamma,1.,0.1,$Par)
end
test(Par)
##
