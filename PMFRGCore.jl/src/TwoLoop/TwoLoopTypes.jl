struct TwoLoop end

struct TwoLoopParams{F,G<:Geometry} <: PMFRGParams
    System::G
    NumericalParams::NumericalParams{F}
    Options::OptionParams
end

getLoopOrder(P::TwoLoopParams) = 2
getPMFRGMethod(P::TwoLoopParams) = TwoLoop()
generateFileName(Par::TwoLoopParams, arg::String = "") = _generateFileName(Par, "_l2" * arg)



## ______________ State Variables shorthand ______________

"""
Convenience struct containing references to arrays for vertices and their derivative. This does not allocate any memory and is performant. 
"""
struct TwoLoopWorkspace{F,Buff,ParType} <: PMFRGWorkspace
    State::StateType{F} #Stores the current state
    Deriv::StateType{F} #Stores the derv

    X::BubbleType{F} #Stores the bubble function X and XTilde
    Y::BubbleType{F} #Stores the bubble function X and XTilde

    Buffer::Buff #Buffer Arrays

    Par::ParType # Params
end

function TwoLoopWorkspace(Deriv::AbstractVector, State::AbstractVector, X, Y, Buffer, Par)
    setZero!(Deriv)
    setZero!(X)
    setZero!(Y)
    return TwoLoopWorkspace(
        StateType(unpackStateVector(State, Par)...),
        StateType(unpackStateVector(Deriv, Par)...),
        X,
        Y,
        Buffer,
        Par,
    )
end

struct BubbleBufferType{T}
    XTa21::Vector{T}

    XTb21::Vector{T}

    XTc21::Vector{T}

    XTd21::Vector{T}
end
BubbleBufferType(type, Npairs) = BubbleBufferType((zeros(type, Npairs) for _ = 1:4)...)

struct BufferTypeTwoLoop{PropsBuff,VertexBuff,XBuff}
    Props::PropsBuff
    Vertex::VertexBuff
    X::XBuff
end
