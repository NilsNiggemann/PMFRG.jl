module StateLib
import PMFRGCore: getVDims, PMFRGParams, StateType

export getF_int,
    getGamma,
    getVa,
    getVb,
    getVc,
    createStateVector,
    unpackStateVector,
    unpackStateVector!,
    getArrayGeometry,
    repackStateVector,
    repackStateVector!


function createStateVector(array_geometry::NamedTuple; floattype)
    zeros(
        floattype,
        array_geometry.NUnique + # f_int
        array_geometry.NUnique * array_geometry.Ngamma + # gamma
        3 * prod(array_geometry.VDims),  # Va, Vb, Vc
    )
end

function getF_int(state::AbstractVector, array_geometry::NamedTuple)
    @view state[_getF_intRange(array_geometry)]
end

function getGamma(state::AbstractVector, array_geometry::NamedTuple)
    v = view(state, _getGammaRange(array_geometry))
    reshape(v, (array_geometry.NUnique, array_geometry.Ngamma))
end

function getVa(state::AbstractVector, array_geometry::NamedTuple)
    v = view(state, _getVaRange(array_geometry))
    reshape(v, array_geometry.VDims)
end

function getVb(state::AbstractVector, array_geometry::NamedTuple)
    v = view(state, _getVbRange(array_geometry))
    reshape(v, array_geometry.VDims)
end

function getVc(state::AbstractVector, array_geometry::NamedTuple)
    v = view(state, _getVcRange(array_geometry))
    reshape(v, array_geometry.VDims)
end


function unpackStateVector(state::AbstractVector, array_geometry::NamedTuple)
    [f(state, array_geometry) for f in [getF_int, getGamma, getVa, getVb, getVc]]
end

function unpackStateVector!(State::StateType, StateArr::AbstractVector)
    array_geometry = getArrayGeometry(State)

    funcs = (_getF_intRange, _getGammaRange, _getVaRange, _getVbRange, _getVcRange)
    parts = (State.f_int, State.γ, State.Γ.a, State.Γ.b, State.Γ.c)

    for (func, part) in zip(funcs, parts)
        reshape(part, :) .= StateArr[func(array_geometry)]
    end

    nothing
end



function getArrayGeometry(Par::PMFRGParams)
    (
        Ngamma = Par.NumericalParams.Ngamma,
        NUnique = Par.System.NUnique,
        VDims = getVDims(Par),
    )
end

function getArrayGeometry(State::StateType)
    (NUnique = size(State.f_int)[1], Ngamma = size(State.γ)[2], VDims = size(State.Γ.a))
end

function unpackStateVector(state::AbstractVector, Par::PMFRGParams)
    unpackStateVector(state, getArrayGeometry(Par))
end

function repackStateVector!(Unrolled::AbstractVector{T}, State::StateType{T}) where {T}
    array_geometry = getArrayGeometry(State)


    funcs = (_getF_intRange, _getGammaRange, _getVaRange, _getVbRange, _getVcRange)
    parts = (State.f_int, State.γ, State.Γ.a, State.Γ.b, State.Γ.c)


    for (func, part) in zip(funcs, parts)
        Unrolled[func(array_geometry)] .= reshape(part, :)
    end

    nothing
end

function repackStateVector!(
    Unrolled::AbstractVector{T},
    r::Base.UnitRange,
    State::StateType{T},
) where {T}
    array_geometry = getArrayGeometry(State)

    funcs = (_getF_intRange, _getGammaRange, _getVaRange, _getVbRange, _getVcRange)
    parts = (State.f_int, State.γ, State.Γ.a, State.Γ.b, State.Γ.c)

    intersections = [intersect(r, f(array_geometry)) for f in funcs]

    intersections_destination =
        [intersection .- (r.start - 1) for intersection in intersections]
    intersections_source = [
        intersection .- (f(array_geometry).start - 1) for
        (intersection, f) in zip(intersections, funcs)
    ]

    for (intersection_destination, part, intersection_source) in
        zip(intersections_destination, parts, intersections_source)
        Unrolled[intersection_destination] .= reshape(part, :)[intersection_source]
    end

    nothing
end



function repackStateVector(State::StateType{T}) where {T}
    array_geometry = getArrayGeometry(State)
    packed = createStateVector(array_geometry; floattype = T)
    repackStateVector!(packed, State)
    packed
end


_getF_intRange(array_geometry::NamedTuple) = 1:array_geometry.NUnique

function _getGammaRange(array_geometry::NamedTuple)
    start = _getF_intRange(array_geometry).stop + 1
    stop = start + array_geometry.NUnique * array_geometry.Ngamma - 1
    start:stop
end

function _getVaRange(array_geometry::NamedTuple)
    start = _getGammaRange(array_geometry).stop + 1
    stop = start + prod(array_geometry.VDims) - 1
    start:stop
end

function _getVbRange(array_geometry::NamedTuple)
    start = _getVaRange(array_geometry).stop + 1
    stop = start + prod(array_geometry.VDims) - 1
    start:stop
end


function _getVcRange(array_geometry::NamedTuple)
    start = _getVbRange(array_geometry).stop + 1
    stop = start + prod(array_geometry.VDims) - 1
    start:stop
end

end
