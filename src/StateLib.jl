module StateLib
import PMFRG: getVDims, PMFRGParams, StateType

export get_f_int, get_gamma, get_Va, get_Vb, get_Vc, CreateState, unpack_state_vector, get_array_geometry, repack!, repack


function CreateState(array_geometry::NamedTuple;floattype)
    zeros(floattype,
        array_geometry.NUnique  + # f_int
        array_geometry.NUnique * array_geometry.Ngamma + # gamma
        3*prod(array_geometry.VDims)  # Va, Vb, Vc
)
end

function get_f_int(state::AbstractVector,array_geometry::NamedTuple)
    @view state[_get_f_int_range(array_geometry)]
end

function get_gamma(state::AbstractVector,array_geometry::NamedTuple)
    v = view(state,_get_gamma_range(array_geometry))
    reshape(v,(array_geometry.NUnique,array_geometry.Ngamma))
end

function get_Va(state::AbstractVector,array_geometry::NamedTuple)
    v = view(state,_get_Va_range(array_geometry))
    reshape(v,array_geometry.VDims)
end

function get_Vb(state::AbstractVector,array_geometry::NamedTuple)
    v = view(state,_get_Vb_range(array_geometry))
    reshape(v,array_geometry.VDims)
end

function get_Vc(state::AbstractVector,array_geometry::NamedTuple)
    v = view(state,_get_Vc_range(array_geometry))
    reshape(v,array_geometry.VDims)
end


function unpack_state_vector(state::AbstractVector,array_geometry::NamedTuple)
    [f(state,array_geometry) for f in [get_f_int,
                                    get_gamma,
                                    get_Va,
                                    get_Vb,
                                    get_Vc]]
end

function get_array_geometry(Par::PMFRGParams)
    (Ngamma = Par.NumericalParams.Ngamma,
     NUnique = Par.System.NUnique,
     VDims = getVDims(Par))
end

function get_array_geometry(State::StateType)
    (NUnique = size(State.f_int)[1],
     Ngamma = size(State.γ)[2],
     VDims = size(State.Γ.a))
end

function unpack_state_vector(state::AbstractVector,Par::PMFRGParams)
    unpack_state_vector(state,get_array_geometry(Par))

end

function repack!(Unrolled::AbstractVector{T},Deriv::StateType{T}) where T
    array_geometry = get_array_geometry(Deriv)
    get_f_int(Unrolled,array_geometry) .= Deriv.f_int
    get_gamma(Unrolled,array_geometry) .= Deriv.γ
    get_Va(Unrolled,array_geometry) .= Deriv.Γ.a
    get_Vb(Unrolled,array_geometry) .= Deriv.Γ.b
    get_Vc(Unrolled,array_geometry) .= Deriv.Γ.c
end

function repack(Deriv::StateType{T}) where T
    array_geometry = get_array_geometry(Deriv)
    packed = CreateState(array_geometry; floattype = T)
    repack!(packed,Deriv)
end


_get_f_int_range(array_geometry::NamedTuple) =  1:array_geometry.NUnique

function _get_gamma_range(array_geometry::NamedTuple)
    start = _get_f_int_range(array_geometry).stop + 1
    stop = start + array_geometry.NUnique * array_geometry.Ngamma - 1
    start:stop
end

function _get_Va_range(array_geometry::NamedTuple)
    start = _get_gamma_range(array_geometry).stop + 1
    stop = start + prod(array_geometry.VDims) - 1
    start:stop
end

function _get_Vb_range(array_geometry::NamedTuple)
    start =  _get_Va_range(array_geometry).stop + 1
    stop = start + prod(array_geometry.VDims) -1
    start:stop
end


function _get_Vc_range(array_geometry::NamedTuple)
    start = _get_Vb_range(array_geometry).stop + 1
    stop = start + prod(array_geometry.VDims) - 1
    start:stop
end

end
