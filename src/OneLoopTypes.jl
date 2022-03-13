struct OneLoop end
struct OptionParams <: AbstractOptions
    usesymmetry::Bool
    MinimalOutput::Bool
end
OptionParams(;usesymmetry::Bool = true,MinimalOutput::Bool = false,kwargs...) = OptionParams(usesymmetry,MinimalOutput)
struct OneLoopParams{F,G <: Geometry} <: PMFRGParams
    System::G
    NumericalParams::NumericalParams{F}
    Options::OptionParams
end

"""Todo: make this error when an unknown kwarg is given!"""
Params(System::Geometry,O::OneLoop=OneLoop();kwargs...) = OneLoopParams(System,NumericalParams(;kwargs...),OptionParams(;kwargs...))


getLoopOrder(P::OneLoopParams) = 1
getLoopMethod(P::OneLoopParams) = OneLoop()

function getAllParamFields(O::Par) where Par #<: PMFRGParams
    fields = Symbol[]
    for f in fieldnames(Par)
        println(getfield(O,f))
        # T = typeof( getfield(O,f))
        # println(T)
        # append!(fields,fieldnames(T))
    end
    return fields
end

## ______________ State Variables shorthand ______________

struct BufferType{PropsBuff,VertexBuff}
    Props::Vector{PropsBuff}
    Vertex::Vector{VertexBuff}
end

"""
Convenience struct containing references to arrays for vertices and their derivative. This does not allocate any memory and is performant. 
"""
struct OneLoopWorkspace{F,PropsBuff,VertexBuff,ParType <: OneLoopParams} <:PMFRGWorkspace
    State::StateType{F} #Stores the current state
    Deriv::StateType{F} #Stores the derv
    
    X::BubbleType{F} #Stores the bubble function X and XTilde

    Buffer::BufferType{PropsBuff,VertexBuff}
 
    Par::ParType # Params

end

function OneLoopWorkspace(Deriv::ArrayPartition,State::ArrayPartition,X,Buffer,Par)
    setZero!(Deriv)
    setZero!(X)
    return OneLoopWorkspace(
        StateType(State.x...),
        StateType(Deriv.x...),
        X,
        Buffer,
        Par
    )
end
