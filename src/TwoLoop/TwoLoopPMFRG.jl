"""
includes Two-loop corrections to PMFRG.
"""

include("TwoLoopTypes.jl")

function AllocateSetup(Par::TwoLoopParams)
    @unpack N,Ngamma = Par.NumericalParams

    @unpack couplings,NUnique,Npairs = Par.System
    println("Two Loop: T= ",Par.NumericalParams.T)
    ##Allocate Memory:
    X = BubbleType(Par)
    Y = BubbleType(Par)
    PropsBuffers = [MMatrix{NUnique,NUnique,double}(undef) for _ in 1:Threads.nthreads()] 
    VertexBuffers = [VertexBufferType(Npairs) for _ in 1:Threads.nthreads()]
    BubbleBuffers = [BubbleBufferType(Npairs) for _ in 1:Threads.nthreads()]
    Buffs = BufferTypeTwoLoop(PropsBuffers,VertexBuffers,BubbleBuffers) 
    return (X,Y,Buffs,Par)
end
include("Buffers.jl")
include("Bubbles.jl")
include("Flowequations.jl")

Params(System::Geometry,O::TwoLoop;kwargs...) = TwoLoopParams(System,NumericalParams(;kwargs...),OptionParams(;kwargs...))

"""
Solves FRG as specified for parameters
"""
SolveFRG(Par::TwoLoopParams;kwargs...) = launchPMFRG!(InitializeState(Par),AllocateSetup(Par),getDeriv2L!; kwargs...)

SolveFRG_Checkpoint(Filename::String,Geometry::SpinFRGLattices.Geometry,Method::TwoLoop,Par=nothing;kwargs...)= launchPMFRG_Checkpoint(Filename,Geometry,AllocateSetup,getDeriv2L!,Par;kwargs...)

generateFileName(Par::TwoLoopParams,Method::TwoLoop,arg::String = "") = generateFileName(Par,"_l2"*arg)