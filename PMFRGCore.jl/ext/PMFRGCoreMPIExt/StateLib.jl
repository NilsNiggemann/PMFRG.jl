function PMFRGCore.repackStateVector!(
    Unrolled::PencilArray{T},
    State::PMFRGCore.StateType{T},
) where {T}
    range = range_local(Unrolled)[1]
    PMFRGCore.repackStateVector!(Unrolled, range, State)
end
