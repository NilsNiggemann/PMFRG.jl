
"""Given set of parameters solve the self-consistent Parquet approximation iteratively.
TODO: Save output to file"""
function SolveParquet(Par::PMFRGCore.ParquetParams, Lam::Real; kwargs...)
    Workspace = PMFRGCore.SetupParquet(Par)
    SolveParquet(Workspace, Lam; kwargs...)
end
function SolveParquet(
    State::PMFRGCore.StateType,
    Par::PMFRGCore.ParquetParams,
    Lam::Real;
    kwargs...,
)
    Workspace = PMFRGCore.SetupParquet(Par)
    writeTo!(Workspace.OldState, State)
    SolveParquet(Workspace, Lam; kwargs...)
end

function SolveParquet(
    Workspace::PMFRGCore.ParquetWorkspace,
    Lam::Real;
    iterator = PMFRGCore.iterateSolution_FP!,
    MainFile = nothing,
    Group = PMFRGCore.DefaultGroup(Workspace.Par),
    CheckpointDirectory = nothing,
    ObsType = PMFRGCore.Observables,
    kwargs...,
)
    Obs = StructArray(ObsType[])
    @time Workspace, Obs = iterator(Workspace, Lam, Obs)
    saveMainOutput(MainFile, Workspace.State, Obs, Lam, Workspace.Par, Group)
    return Workspace, Obs
end
