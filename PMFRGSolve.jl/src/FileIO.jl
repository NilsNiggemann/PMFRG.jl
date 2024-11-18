using PMFRGCore, SpinFRGLattices, HDF5, DiffEqCallbacks

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
        saveMainOutput(MainFile, sol, saved_values_full, Par, Group)
    end
    return sol, saved_values_full
end

SolveFRG_Checkpoint(
    Filename::String,
    GeometryGenerator::Function,
    Par = nothing;
    kwargs...,
) = SolveFRG_Checkpoint(Filename, readGeometry(Filename, GeometryGenerator), Par; kwargs...)

function saveCurrentState(
    DirPath::String,
    State::AbstractArray,
    saved_Values::DiffEqCallbacks.SavedValues,
    Lam::Real,
    Par::PMFRGParams,
)
    Filename = joinpath(DirPath, "CurrentState.h5")
    saveState(Filename, State, Lam, Par, "w")
    saveParams(Filename, Par)
    saveObs(Filename, saved_Values, "Observables")
    Filename
end
saveCurrentState(::Nothing, args...) = nothing



saveMainOutput(
    Filename::String,
    Solution::ODESolution,
    saved_values::DiffEqCallbacks.SavedValues,
    Par::PMFRGParams,
    Group::String,
) = saveMainOutput(Filename, Solution.u[end], saved_values, saved_values.t[end], Par, Group)

function saveMainOutput(
    Filename::String,
    State,
    saved_values,
    Lambda::Real,
    Par::PMFRGParams,
    Group::String,
)
    function save(name)
        saveMainOutput(name, saved_values, Group)
        saveExtraFields(name, State, Lambda, Par, Group)
    end
    try
        save(Filename)
    catch e
        newName = UniqueFileName(Filename)
        @warn "Writing to $Filename errored with exception $(string(e))! Writing to $newName instead."
        save(newName)
    end
    # saveParams(Filename,Par,Group)
end

saveMainOutput(
    Filename::String,
    Solution::ODESolution,
    saved_values,
    Par::PMFRGParams,
    Group::Nothing,
) = saveMainOutput(Filename, Solution, saved_values, Par, string(Par.NumericalParams.T))

function saveMainOutput(Filename::String, saved_values, Group::String)
    mkpath(dirname(Filename))
    println("Saving Main output to ", abspath(Filename))
    saveObs(Filename, saved_values, Group)
end
saveMainOutput(::Nothing, args...) = nothing


function setCheckpoint(
    Directory::String,
    State,
    saved_values::DiffEqCallbacks.SavedValues,
    Lam,
    Par,
    checkPointList,
)
    saveCurrentState(Directory, State, saved_values, Lam, Par)
    if !isempty(checkPointList)
        if Lam < last(checkPointList)
            Checkpoint = pop!(checkPointList)
            CheckpointFile = UniqueFileName("$Directory/$(strd(Lam)).h5")
            println("\nsaving Checkpoint Lam ≤ $Checkpoint at ", Lam)
            println("in file ", CheckpointFile)
            println("")
            mv(joinpath(Directory, "CurrentState.h5"), CheckpointFile)
        end
    end
end

function setCheckpoint(Directory::Nothing, State, saved_values, Lam, Par, checkPointList)
    return
end

"""Saves Observables"""
function saveObs(
    Filename::String,
    saved_values::DiffEqCallbacks.SavedValues,
    Group::String = "",
)
    ObsArr = StructArray(saved_values.saveval)
    saveObs(Filename, ObsArr, Group)
    h5write(Filename, PMFRGCore.joinGroup(Group, "Lambda"), saved_values.t)
end

function saveObs(Filename::String, Obs::StructArray{ObsType}, Group::String) where {ObsType}

    function convertToArray(
        VecOfArray::AbstractVector{VT},
    ) where {N,VT<:AbstractArray{T,N} where {T}}
        cat(VecOfArray..., dims = N + 1)
    end
    AbsoluteFilename = abspath(Filename)

    Fields = fieldnames(ObsType)
    for F in Fields
        arr = convertToArray(getproperty(Obs, F))
        h5write(AbsoluteFilename, PMFRGCore.joinGroup(Group, string(F)), arr)
    end
end

function readObservables(Filename::String, ObsType = PMFRGCore.Observables)
    t = h5read(Filename, "Observables/Lambda")
    Fields = fieldnames(ObsType)
    obsTuple = Tuple(h5read(Filename, "Observables/$f") for f in Fields)
    saveval = ObsType[]
    saved_values = SavedValues(eltype(t), ObsType)
    # return obsTuple
    for i in eachindex(t)
        Currentval(x) = selectdim(x, length(size(x)), i) |> Array
        # return Currentval.(obsTuple)
        push!(saveval, ObsType(Currentval.(obsTuple)...))
    end
    append!(saved_values.t, t)
    append!(saved_values.saveval, saveval)
    return saved_values
end

function saveExtraFields(
    Filename::String,
    State,
    Lambda::Real,
    Par::PMFRGParams,
    Group::String,
)
    (; T, N) = Par.NumericalParams
    (; Name, NUnique, Npairs, NLen) = Par.System
    Chi_nu = getChi(State, Lambda, Par, N)
    h5write(Filename, "$Group/Name", Name)
    h5write(Filename, "$Group/Npairs", Npairs)
    h5write(Filename, "$Group/T", T)
    h5write(Filename, "$Group/N", N)
    h5write(Filename, "$Group/NUnique", NUnique)
    h5write(Filename, "$Group/NLen", NLen)
    h5write(Filename, "$Group/Chi_nu", Chi_nu)
end
