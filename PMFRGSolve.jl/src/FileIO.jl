using PMFRGCore, SpinFRGLattices, HDF5

# Todo: provide SpinFRGLattices.getGeometryGenerator that takes Name string and returns correct method
SolveFRG_Checkpoint(
    Filename::String,
    Geometry::SpinFRGLattices.Geometry,
    Par = nothing;
    kwargs...,
) = launchPMFRG_Checkpoint(
    Filename,
    Geometry,
    PMFRGCore.AllocateSetup,
    PMFRGCore.getDeriv!,
    Par;
    kwargs...,
)


function launchPMFRG_Checkpoint(
    Filename::String,
    Geometry::SpinFRGLattices.Geometry,
    AllocatorFunction::Function,
    Derivative::Function,
    Par = nothing;
    MainFile = nothing,
    Group = nothing,
    Params = (),
    ObservableType = PMFRGCore.Observables,
    kwargs...,
)
    State = readState(Filename)
    Old_Lam_max = h5read(Filename, "Params/Lam_max")
    Par = PMFRGCore.getFileParams(Filename, Geometry, Par; Params...)
    saved_values_full = readObservables(Filename, ObservableType)
    setup = AllocatorFunction(Par)
    CheckPointfolder = dirname(Filename)
    FilePath = dirname(CheckPointfolder)
    ObsSaveat = PMFRGCore.getLambdaMesh(nothing, Par.NumericalParams.Lam_min, Old_Lam_max)
    filter!(x -> x < Par.NumericalParams.Lam_max, ObsSaveat)
    sol, saved_values = launchPMFRG!(
        State,
        setup,
        Derivative;
        CheckpointDirectory = FilePath,
        ObsSaveat = ObsSaveat,
        kwargs...,
        MainFile = nothing,
    ) #launch PMFRG but do not save output yet
    append!(saved_values_full.t, saved_values.t)
    append!(saved_values_full.saveval, saved_values.saveval)
    if MainFile !== nothing
        PMFRGCore.saveMainOutput(MainFile, sol, saved_values_full, Par, Group)
    end
    return sol, saved_values_full
end

SolveFRG_Checkpoint(
    Filename::String,
    GeometryGenerator::Function,
    Par = nothing;
    kwargs...,
) = SolveFRG_Checkpoint(Filename, readGeometry(Filename, GeometryGenerator), Par; kwargs...)
