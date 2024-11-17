using PencilArrays

"""
Collect the state information from different MPI ranks
and returns a new global array.
"""
function gatherState(LocalState::PencilArray, Par::PMFRGCore.PMFRGParams)

    GlobalState = createGlobalStateVector(Par)
    gatherState!(GlobalState, LocalState)
    GlobalState

end

"""
Collect the state information from different MPI ranks
into an existing GlobalState vector.
"""
function gatherState!(GlobalState::AbstractVector{T}, LocalState::PencilArray{T}) where {T}

    GlobalStateBuff = VBuffer(
        GlobalState,
        let
            counts = map(
                r -> length(range_remote(LocalState, r + 1)[1]),
                LocalState.pencil.topology.ranks,
            )
        end,
    )

    gatherState!(GlobalStateBuff, LocalState)

end

"""
Collect the state information from different MPI ranks
into an existing VBuffer objest, that references a GlobalState.
"""
function gatherState!(GlobalStateBuff::VBuffer, LocalState::PencilArray)

    GlobalStateVector = GlobalStateBuff.data
    GlobalStateVector[range_local(LocalState)[1]] = LocalState
    MPI.Allgatherv!(GlobalStateBuff, MPI.COMM_WORLD)
end


function PMFRGCore.repackStateVector!(
    Unrolled::PencilArray{T},
    State::PMFRGCore.StateType{T},
) where {T}
    range = range_local(Unrolled)[1]
    PMFRGCore.repackStateVector!(Unrolled, range, State)
end

function PMFRGCore.unpackStateVector(LocalState::PencilArray, par::PMFRGCore.PMFRGParams)
    GlobalState = gatherState(LocalState, par)
    PMFRGCore.unpackStateVector(GlobalState, par)
end
