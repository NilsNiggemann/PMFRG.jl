"""
This file contains functions used for reading and writing FRG data to disc.
Stored data can be used to re-launch an incomplete FRG calculation.
"""
joinGroup(args...) = join(args, "/")

"""Saves Vertices to a compressed HDF5 file in a Group "Lam"."""
function saveState(Filename::String, State::ArrayPartition, Lam, mode = "cw")
    Vertices = State.x
    Names = "fint", "gamma", "Va", "Vb", "Vc"
    # Filename = string(DirName,"/$(string(round(Lam,digits =3))).h5")
    try
        h5open(Filename, mode) do f
            for (Name, V) in zip(Names, Vertices)
                f["$Name", blosc = 9] = V
            end
        end
        h5write(Filename, "Lam", Lam)
    catch e
        @warn "Saving state was unsuccessfull with exception:\n $e "
    end
end

"""Reads Vertices from file"""
function readState(Filename::String)
    Names = "fint", "gamma", "Va", "Vb", "Vc"
    State = ArrayPartition((h5read(Filename, "$N") for N in Names)...)
    return State
end

function readLam(Filename::String)
    Lam = h5read(Filename, "Lam")
    return Lam
end

function readObservables(Filename::String, ObsType = Observables)
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

function readGeometry(Filename::String, GeometryGenerator::Function; Group = "Geometry")
    NLen = h5read(Filename, "$Group/NLen")
    couplings = h5read(Filename, "$Group/couplings")
    Npairs = h5read(Filename, "$Group/Npairs")
    Name = h5read(Filename, "$Group/Name")
    System = GeometryGenerator(NLen)
    @assert Npairs == System.Npairs "Number of unique pairs does not match geometry in $Filename"
    @assert Name == System.Name "Name does not match geometry in $Filename"
    System.couplings .= couplings
    return System
end

EssentialParamFields() =
    (:T, :N, :Ngamma, :accuracy, :Lam_min, :Lam_max, :ex_freq, :lenIntw, :lenIntw_acc)

"""Saves important information about computation parameters so that they can be reconstructed"""
function saveNumericalParams(Filename, Par::PMFRGParams, Group = "")
    Fields = EssentialParamFields()
    for F in Fields
        h5write(Filename, joinGroup(Group, "Params/$F"), getfield(Par.NumericalParams, F))
    end
end

"""Saves important information about Geometry parameters so that they can be reconstructed"""
function saveGeometryParams(Filename, Par::PMFRGParams, Group = "")
    h5write(Filename, joinGroup(Group, "Geometry/Name"), Par.System.Name)
    h5write(Filename, joinGroup(Group, "Geometry/NLen"), Par.System.NLen)
    h5write(Filename, joinGroup(Group, "Geometry/couplings"), Par.System.couplings)
    h5write(Filename, joinGroup(Group, "Geometry/Npairs"), Par.System.Npairs)
end

function saveMethodParams(Filename, Par::PMFRGParams, Group = "")
    h5write(Filename, joinGroup(Group, "Params/looporder"), getLoopOrder(Par))
end
readLoopOrder(Filename, Group = "") = h5read(Filename, string(Group, "/Params/looporder"))


"""Saves important information about parameters so that they can be reconstructed
"""
function saveParams(Filename, Par::PMFRGParams, Group = "")
    saveNumericalParams(Filename, Par, Group)
    saveGeometryParams(Filename, Par, Group)
    saveMethodParams(Filename, Par, Group)
end

function readParams(Filename::String, Geometry::SpinFRGLattices.Geometry; modifyParams...)
    Fields = EssentialParamFields()
    Kwargs = Dict((F => h5read(Filename, "Params/$F") for F in Fields)...)
    Method = getPMFRGMethod(readLoopOrder(Filename))
    Par = Params(Geometry, Method; Kwargs..., modifyParams...)
    return Par
end

readParams(Filename::String, GeometryGenerator::Function; modifyParams...) =
    readParams(Filename, readGeometry(Filename, GeometryGenerator); modifyParams...)

function modifyParams(Par; modifyParams...)
    NumParKwargs = Dict(
        (F => getfield(Par.NumericalParams, F) for F in fieldnames(NumericalParams))...,
    )
    Par = Params(Par.System, getPMFRGMethod(Par); NumParKwargs..., modifyParams...)
end

function setupDirectory(DirPath, Par; overwrite = false)
    DirPath = generateName_verbose(DirPath, Par)
    overwrite || (DirPath = UniqueDirName(DirPath))
    println("Checkpoints saved at $(abspath(DirPath))")
    # CheckPath = joinpath(DirPath,"Checkpoints")
    mkpath(DirPath)
    return DirPath
end

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

"""Rename CurrentState to FinalState as indicator that Job is finished"""
function SetCompletionCheckmark(DirPath::String)
    CurrState = joinpath(DirPath, "CurrentState.h5")
    FinState = UniqueFileName(joinpath(DirPath, "FinalState.h5"))
    if ispath(CurrState)
        mv(CurrState, FinState)
    end
end
SetCompletionCheckmark(DirPath::Nothing) = nothing

file_extension_pos(file::String) = findlast('.', file)
function file_extension(file::String)
    pos = file_extension_pos(file)
    pos === nothing && return nothing
    return file[pos+1:end]
end

function getVersionNumber(Path)
    it = findlast("(v_", Path)
    if it === nothing
        return 0
    else
        it_1 = last(it)
        brackindex = findlast(')', Path)
        try
            return parse(Int, Path[it_1+1:brackindex-1])

        catch e
            println(Path, Path[it_1+1:brackindex-1])
            throw(e)
        end
    end
end

versionString(versNum) = "(v_$versNum)"

function UniqueFileName(Path)
    ending = "." * file_extension(Path)
    while ispath(Path)
        versNum = getVersionNumber(Path)
        if versNum == 0
            Path = replace(Path, ending => versionString(versNum + 1) * ending)
        end
        Path = replace(Path, versionString(versNum) => versionString(versNum + 1))
    end
    return Path
end

function UniqueDirName(Path)
    while ispath(Path)
        versNum = getVersionNumber(Path)
        if versNum == 0
            Path = Path * versionString(versNum + 1)
        end
        Path = replace(Path, versionString(versNum) => versionString(versNum + 1))
    end
    return Path
end

function generateName_verbose(Directory::String, Par::PMFRGParams)
    (; T, N) = Par.NumericalParams
    Name = Par.System.Name
    Name = "$(Name)_N=$(N)_T=$(T)"
    Name = joinpath(Directory, Name)
    return Name
end

generateUniqueName(Directory::String, Par::PMFRGParams) =
    UniqueDirName(generateName_verbose(Directory, Par))


function _generateFileName(Par::PMFRGParams, arg::String = ""; kwargs...)
    Name = Par.System.Name
    N = Par.NumericalParams.N

    FName =
        "PMFRG_$(Name)_N=$(N)$arg" * join("_$(k)=$(strd(v))" for (k, v) in kwargs) * ".h5"

    return FName
end



generateFileName(Par::PMFRGParams, arg::String = ""; kwargs...) =
    _generateFileName(Par, arg; kwargs...)
generateFileName(Par::OneLoopParams, arg::String = ""; kwargs...) =
    _generateFileName(Par, "_l1" * arg; kwargs...)

function generateMainFile(
    Directory::String,
    Par::PMFRGParams,
    arg::String = "";
    savekeywords = true,
    unique = true,
    kwargs...,
)
    FName = joinpath(Directory, generateFileName(Par, arg; kwargs...))

    unique && (FName = UniqueFileName(FName))
    h5open(FName, "cw") do f
        if savekeywords
            for (k, v) in kwargs
                f[string(k)] = v
            end
        end
    end
    return FName
end
generateMainFile(Par::PMFRGParams, arg::String = ""; kwargs...) =
    generateMainFile(".", Par, arg; kwargs...)

function ParamsCompatible(Par1, Par2)
    Fields = (:Npairs, :N, :Ngamma, :Lam_max, :System)
    for f in Fields
        @assert getproperty(Par1, f) == getproperty(Par2, f) "$f not compatible with parameters used in Checkpoint"
    end

    return true
end

function getFileParams(Filename, Geometry, Par::PMFRGParams)
    Lam = readLam(Filename)
    Par_file = readParams(Filename, Geometry; Lam_max = Lam)
    Par = modifyParams(Par; Lam_max = Lam)
    ParamsCompatible(Par, Par_file)
    return Par
end

function getFileParams(Filename, Geometry, Par::Nothing; kwargs...)
    Lam = readLam(Filename)
    return readParams(Filename, Geometry; Lam_max = Lam, kwargs...)
end
# Todo: provide SpinFRGLattices.getGeometryGenerator that takes Name string and returns correct method
SolveFRG_Checkpoint(
    Filename::String,
    Geometry::SpinFRGLattices.Geometry,
    Par = nothing;
    kwargs...,
) = launchPMFRG_Checkpoint(Filename, Geometry, AllocateSetup, getDeriv!, Par; kwargs...)


function launchPMFRG_Checkpoint(
    Filename::String,
    Geometry::SpinFRGLattices.Geometry,
    AllocatorFunction::Function,
    Derivative::Function,
    Par = nothing;
    MainFile = nothing,
    Group = nothing,
    Params = (),
    ObservableType = Observables,
    kwargs...,
)
    State = readState(Filename)
    Old_Lam_max = h5read(Filename, "Params/Lam_max")
    Par = getFileParams(Filename, Geometry, Par; Params...)
    saved_values_full = readObservables(Filename, ObservableType)
    setup = AllocatorFunction(Par)
    CheckPointfolder = dirname(Filename)
    FilePath = dirname(CheckPointfolder)
    ObsSaveat = getLambdaMesh(nothing, Par.NumericalParams.Lam_min, Old_Lam_max)
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

function saveObs(Filename::String, Obs::StructArray{ObsType}, Group::String) where {ObsType}
    Fields = fieldnames(ObsType)
    for F in Fields
        arr = convertToArray(getproperty(Obs, F))
        h5write(Filename, joinGroup(Group, string(F)), arr)
    end
end


function convertToArray(
    VecOfArray::AbstractVector{VT},
) where {N,VT<:AbstractArray{T,N} where {T}}
    cat(VecOfArray..., dims = N + 1)
end

function saveMainOutput(Filename::String, saved_values, Group::String)
    mkpath(dirname(Filename))
    println("Saving Main output to ", abspath(Filename))
    saveObs(Filename, saved_values, Group)
end
saveMainOutput(::Nothing, args...) = nothing


function setCheckpoint(Directory::String, State, saved_values, Lam, Par, checkPointList)
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
    saved_values::DiffEqCallbacks.SavedValues,
    Par::PMFRGParams,
    Group::Nothing,
) = saveMainOutput(Filename, Solution, saved_values, Par, string(Par.NumericalParams.T))

function getFilesFromSubDirs(Folder::String)
    allpaths = collect(walkdir(Folder))
    Files = String[]
    for p in allpaths
        for filename in p[end]
            pathAndName = joinpath(p[begin], filename)
            push!(Files, pathAndName)
        end
    end
    return Files
end

function getUnfinishedJobs(Folder::String)
    allFiles = filter!(x -> occursin("CurrentState.h5", x), getFilesFromSubDirs(Folder))
end
