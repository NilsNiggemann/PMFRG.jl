@proto struct OneLoopParams{F} <: PMFRGParams
    System::Geometry
    NumericalParams::NumericalParams{F}
    Options::OptionParams
end

getLoopOrder(P::OneLoopParams) = 1


function OneLoopParams(System;NumericalParams,Options,kwargs..) 
    OneLoopParams(System,NumericalParams(kwargs...),Options(kwargs...))
end

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
"""Struct containing the observables that are saved at every step"""
struct Observables{T}
    Chi::Vector{T}
    gamma::Matrix{T}
    f_int::Vector{T}
    MaxVa::Vector{T}
    MaxVb::Vector{T}
    MaxVc::Vector{T}
end

## ______________ State Variables shorthand ______________

"""
Convenience struct containing references to arrays for vertices and their derivative. This does not allocate any memory and is performant. 
"""
struct OneLoopWorkspace{F,PropsBuff,VertexBuff}
    State::StateType{F} #Stores the current state
    Deriv::StateType{F} #Stores the derv
    
    X::BubbleType{F} #Stores the bubble funtion X and XTilde

    Par::OneLoopParams{F} # Params

    Buffer::BufferType{PropsBuff,VertexBuff}

end

function OneLoopWorkspace(Deriv::ArrayPartition,State::ArrayPartition,X,Par)
    setZero!(Deriv)
    setZero!(X)
    return OneLoopWorkspace(
        VertexType(State.x...),
        VertexType(Deriv.x...),
        X,
        Par
    )
end
# given a Matsubara (!) integer, return the corresponding Matsubara frequency
function get_w(nw,T)
    return (2*nw+1)*pi*T
end

function absmax(x)
    MaxAndPos = findmax(x)
    MinAndPos = findmin(x)
    absvals = (abs(MaxAndPos[1]),abs(MinAndPos[1]))
    return (MaxAndPos,MinAndPos)[argmax(absvals)]
end
strd(x) = string(round(x,digits=3))
