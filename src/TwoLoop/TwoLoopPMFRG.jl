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
    floattype = _getFloatType(Par)
    PropsBuffers = [MMatrix{NUnique,NUnique,floattype}(undef) for _ in 1:Threads.nthreads()] 
    VertexBuffers = [VertexBufferType(floattype,Npairs) for _ in 1:Threads.nthreads()]
    BubbleBuffers = [BubbleBufferType(floattype,Npairs) for _ in 1:Threads.nthreads()]
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
SolveFRG(Par::TwoLoopParams;kwargs...) = launchPMFRG!(InitializeState(Par),AllocateSetup(Par),getDeriv!; kwargs...)

generateFileName(Par::TwoLoopParams,Method::TwoLoop,arg::String = "") = _generateFileName(Par,"_l2"*arg)