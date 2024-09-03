struct OneLoop end
struct OptionParams <: AbstractOptions
    usesymmetry::Bool
    MinimalOutput::Bool
end

OptionParams(; usesymmetry::Bool = true, MinimalOutput::Bool = false, kwargs...) =
    OptionParams(usesymmetry, MinimalOutput)

abstract type AbstractOneLoopParams <: PMFRGParams end


struct OneLoopParams{F<:AbstractFloat,G<:Geometry} <: AbstractOneLoopParams
    System::G
    NumericalParams::NumericalParams{F}
    Options::OptionParams
end


"""Construct a set of parameters used for an FRG calculation. 'System' refers to the geometry that is used for the calculation (see SpinFRGLattices).
The second argument is the method, i.e. OneLoop() (default value), TwoLoop(), Parquet(). For possible optional keyword arguments see the docs of 'NumericalParams'.
"""
Params(System::Geometry, ::OneLoop = OneLoop(); kwargs...) =
    OneLoopParams(System, NumericalParams(; kwargs...), OptionParams(; kwargs...)) # Todo: make this error when an unknown kwarg is given!

getLoopOrder(P::AbstractOneLoopParams) = 1
getPMFRGMethod(P::AbstractOneLoopParams) = OneLoop()

function getAllParamFields(O::Par) where {Par} #<: PMFRGParams
    fields = Symbol[]
    for f in fieldnames(Par)
        println(getfield(O, f))
        # T = typeof( getfield(O,f))
        # println(T)
        # append!(fields,fieldnames(T))
    end
    return fields
end

## ______________ State Variables shorthand ______________

struct BufferType{PropsBuff,VertexBuff}
    Props::PropsBuff
    Vertex::VertexBuff
end

"""
Convenience struct containing references to arrays for vertices and their derivative. This does not allocate any memory and is performant. 
"""
struct OneLoopWorkspace{F,PropsBuff,VertexBuff,ParType<:AbstractOneLoopParams} <:
       PMFRGWorkspace
    State::StateType{F} #Stores the current state
    Deriv::StateType{F} #Stores the derivative

    X::BubbleType{F} #Stores the bubble function X and XTilde

    Buffer::BufferType{PropsBuff,VertexBuff}

    Par::ParType # Params

end

function OneLoopWorkspace(Deriv::ArrayPartition, State::ArrayPartition, X, Buffer, Par)
    setZero!(Deriv)
    setZero!(X)
    return OneLoopWorkspace(StateType(State.x...), StateType(Deriv.x...), X, Buffer, Par)
end
