"""
includes Two-loop corrections to PMFRG.
"""

include("TwoLoopTypes.jl")

function AllocateSetup(
    Par::TwoLoopParams,
    ParallelizationScheme::AbstractParallelizationScheme = MultiThreaded(),
)

    (; NUnique, Npairs) = Par.System
    Par.Options.MinimalOutput || println("Two Loop: T= ", Par.NumericalParams.T)
    ##Allocate Memory:
    X = BubbleType(Par)
    Y = BubbleType(Par)
    floattype = _getFloatType(Par)
    PropsBuffers =
        getChannel([MMatrix{NUnique,NUnique,floattype}(undef) for _ = 1:Threads.nthreads()])
    VertexBuffers =
        getChannel([VertexBufferType(floattype, Npairs) for _ = 1:Threads.nthreads()])
    BubbleBuffers =
        getChannel([BubbleBufferType(floattype, Npairs) for _ = 1:Threads.nthreads()])
    Buffs = BufferTypeTwoLoop(PropsBuffers, VertexBuffers, BubbleBuffers)
    return (; X, Y, Buffs, Par, ParallelizationScheme)
end
include("Buffers.jl")
include("Bubbles.jl")
include("Flowequations.jl")

Params(System::Geometry, O::TwoLoop; kwargs...) =
    TwoLoopParams(System, NumericalParams(; kwargs...), OptionParams(; kwargs...))

"""
Solves FRG as specified for parameters
"""
SolveFRG(Par::TwoLoopParams; kwargs...) =
    launchPMFRG!(InitializeState(Par), AllocateSetup(Par), getDerivTwoLoop!; kwargs...)

generateFileName(Par::TwoLoopParams, Method::TwoLoop, arg::String = "") =
    _generateFileName(Par, "_l2" * arg)
