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

    Par = Params(System = System,N=24,T=1.2,MinimalOutput=false,usesymmetry = true,accuracy= 1E-4,Ngamma = 24,Lam_min = 10)
    Solution,saved_values = SolveFRG(Par,TwoLoop(),MainFile = "/storage/niggeni/test/main.h5",CheckpointDirectory = "/storage/niggeni/test/",method = BS3(),VertexCheckpoints = [],overwrite_Checkpoints=true)
    # Solution,saved_values = SolveFRG_Checkpoint("/storage/niggeni/test/Cubic_NLen=4_N=24_T=1.200000000000000/CurrentState.h5",getCubic,OneLoop(),MainFile="/storage/niggeni/test/main.h5",method = BS3(),VertexCheckpoints = [],overwrite_Checkpoints = true)
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
using PMFRG,SpinFRGLattices,SpinFRGLattices.Pyrochlore,BenchmarkTools

function CheckKernels(Par::Params;Lam=1.,MaxLengthInSecs=60*1)
    State= PMFRG.InitializeState(Par)
    rand!(State)
    setup = PMFRG.AllocateSetup(Par)
    Deriv = PMFRG.similar(State)

    X,XTilde,PropsBuffers,VertexBuffers,Par = setup #use pre-allocated X and XTilde to reduce garbage collector time
    N = Par.N
    Workspace = PMFRG.Workspace_Struct(Deriv,State,X,XTilde)
    BubbleProp = PropsBuffers[Threads.threadid()] # get pre-allocated thread-safe buffers
    Buffer = VertexBuffers[Threads.threadid()]
    sprop = rand!(BubbleProp)
    display(@benchmark PMFRG.addXTilde!($Workspace,1,1,2,1,$sprop,$Par))
    display(@benchmark PMFRG.addX!($Workspace,1,1,2,1,$sprop,$Par,$Buffer))
end


System = getPyrochlore(5);
Par = Params(System = System,N=24);
State= PMFRG.InitializeState(Par);
rand!(State);
setup = PMFRG.AllocateSetup(Par);
Deriv = PMFRG.similar(State);

X,XTilde,PropsBuffers,VertexBuffers,Par = setup #use pre-allocated X and XTilde to reduce garbage collector time
N = Par.N;
Workspace = PMFRG.Workspace_Struct(Deriv,State,X,XTilde);
BubbleProp = PropsBuffers[Threads.threadid()] # get pre-allocated ;thread-safe buffers
Buffer = VertexBuffers[Threads.threadid()];
sprop = rand!(BubbleProp);

@code_warntype PMFRG.getVertexDeriv!(Workspace,1.,Par,PropsBuffers,VertexBuffers)


##
using PMFRG,SpinFRGLattices
ParquetLambda = 0.
Par = Params(System = getPolymer(2),N=24,T=0.5,accuracy = 1e-5,usesymmetry = true,Lam_min = ParquetLambda)
Sol,Obs = SolveParquet(Par,ParquetLambda,maxiterBubble = 100,maxitergamma = 2000)

##
SolP,ObsPt = SolveFRG(Par)
ObsP = PMFRG.StructArray(ObsPt.saveval)
Obst = ObsPt.t
LamIndex = findfirst(x-> x <= ParquetLambda, Obst)
##
chi1_ex(B) = (exp(B) - 1 + B )/(2*(exp(B)+3))
chi2_ex(B) = -(exp(B) - 1 - B )/(2*(exp(B)+3))
using Plots
numiter = length(Obs)
# plot(PMFRG.convertToArray(Obs.gamma)[1,1:5,:]')
pls = []
plot(PMFRG.convertToArray(Obs.gamma)[1,1:end,2],label = "i=1",ylabel = "γ",xlabel = "n")
plot!(PMFRG.convertToArray(Obs.gamma)[1,1:end,div(numiter,2)],label = "i=$(div(numiter,2))")
push!(pls, plot!(PMFRG.convertToArray(Obs.gamma)[1,1:end,end],label = "i=$(numiter)",color = :black))
plot(xlabel = "iterations",ylabel = "χ")
hline!([chi1_ex(1/Par.T),chi2_ex(1/Par.T)],labels = ["exact" ""],color = :grey, linestyle = :dash)
hline!([ObsP.Chi[LamIndex]'],labels = ["PMFRG" ""],color = :red, linestyle = :dash)
push!(pls, plot!(PMFRG.convertToArray(Obs.Chi)',color = :black,labels = ["Parquet" ""]))

push!(pls, plot(PMFRG.convertToArray(Obs.MaxVc)',ylabel = "max(Γc)",xlabel = "iterations"))

plot(pls...,size = (800,800),title = "T=$(PMFRG.strd(Par.T))")

##
using PMFRG,SpinFRGLattices
ParquetLambda = 0.
chis = Vector{Float64}[]
chisbegin = Vector{Float64}[]
Trange = 0.8:0.1:2.0
for T in Trange
    Par = Params(System = getPolymer(2),N=32,T=T,accuracy = 1e-5,usesymmetry = true,Lam_min = ParquetLambda)
    Sol,Obs = SolveParquet(Par,ParquetLambda,maxiterBubble = 70,maxitergamma = 2000)
    push!(chis,Obs.Chi[end])
    push!(chisbegin,Obs.Chi[begin])
end

##
using LaTeXStrings
bareSt = PMFRG.MultiLoopPMFRG.StateType(Par)
PMFRG.MultiLoopPMFRG.setToBareVertex!(bareSt.Γ,Par)

TrangeDense = 0.0:0.01:maximum(Trange)
TrangeDense2 = 0.5:0.01:maximum(Trange)
chisBare = Vector{Float64}[]
for T in TrangeDense2
    Par = Params(System = getPolymer(2),N=32,T=T,accuracy = 1e-5,usesymmetry = true,Lam_min = ParquetLambda)

    push!(chisBare, PMFRG.MultiLoopPMFRG.getChi(bareSt,Par.Lam_min,Par))
end
plot(xlabel = L"T",ylabel = L"\chi_{ij}")
plot!(TrangeDense,chi1_ex.(1 ./TrangeDense),color = :grey, label = "exact")
plot!(TrangeDense,chi2_ex.(1 ./TrangeDense),color = :grey, label = "")
plot!(TrangeDense2,PMFRG.convertToArray(chisBare)',color = :black, labels = [L"\Gamma = \Gamma^0" ""])
scatter!(Trange,PMFRG.convertToArray(chis)',color = :red,labels = ["Parquet" ""])
# scatter!(Trange,PMFRG.convertToArray(chisbegin)',color = :blue,labels = ["init" ""])

##
@with_kw struct test1
    t::Int
end
test1(;t=1,kwargs...) = test1(t)
@with_kw struct test2
    u::Int
end
test2(;u=2,kwargs...) = test2(u)
struct test3
    a::test1
    b::test2
end
test3(;kwargs...) = test3(test1(;kwargs...),test2(;kwargs...))
a = test3(t=1,u=3)