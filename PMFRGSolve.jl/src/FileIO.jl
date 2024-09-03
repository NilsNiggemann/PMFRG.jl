import PMFRG: PMFRGParams

saveMainOutput(
    Filename::String,
    Solution::ODESolution,
    saved_values::DiffEqCallbacks.SavedValues,
    Par::PMFRGParams,
    Group::String,
) = saveMainOutput(Filename, Solution.u[end], saved_values, saved_values.t[end], Par, Group)

saveMainOutput(
    Filename::String,
    Solution::ODESolution,
    saved_values::DiffEqCallbacks.SavedValues,
    Par::PMFRGParams,
    Group::Nothing,
) = saveMainOutput(Filename, Solution, saved_values, Par, string(Par.NumericalParams.T))

function saveCurrentState(
    DirPath::String,
    State::AbstractArray,
    saved_Values::DiffEqCallbacks.SavedValues,
    Lam::Real,
    Par::PMFRGParams,
)
    Filename = joinpath(DirPath, "CurrentState.h5")
    saveState(Filename, State, Lam, "w")
    saveParams(Filename, Par)
    saveObs(Filename, saved_Values, "Observables")
    Filename
end

saveCurrentState(
    DirPath::Nothing,
    State::AbstractArray,
    saved_Values::DiffEqCallbacks.SavedValues,
    Lam::Real,
    Par::PMFRGParams,
) = nothing

"""Saves Observables"""
function saveObs(
    Filename::String,
    saved_values::DiffEqCallbacks.SavedValues,
    Group::String = "",
)
    ObsArr = StructArray(saved_values.saveval)
    saveObs(Filename, ObsArr, Group)
    h5write(Filename, joinGroup(Group, "Lambda"), saved_values.t)
end
