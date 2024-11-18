"""
This file contains functions used for reading and writing FRG data to disc.
Stored data can be used to re-launch an incomplete FRG calculation.
"""
joinGroup(args...) = join(args, "/")

"""Saves Vertices to a compressed HDF5 file in a Group "Lam"."""
function saveState(
    Filename::String,
    State::AbstractVector,
    Lam,
    Par::PMFRGParams,
    mode = "cw",
)
    Vertices = unpackStateVector(State, Par)
    Names = "fint", "gamma", "Va", "Vb", "Vc"
    # Filename = string(DirName,"/$(string(round(Lam,digits =3))).h5")
    try
        h5open(Filename, mode) do f
            for (Name, V) in zip(Names, Vertices)
                f["$Name", blosc = 9] = Array(V)
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

    f_int, γ, Γa, Γb, Γc = [h5read(Filename, "$N") for N in Names]
    State = StateType(f_int, γ, VertexType(Γa, Γb, Γc))
    StateArray = repackStateVector(State)
    return StateArray
end

function readLam(Filename::String)
    Lam = h5read(Filename, "Lam")
    return Lam
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
    absFilename = abspath(Filename)
    for F in Fields
        h5write(
            absFilename,
            joinGroup(Group, "Params/$F"),
            getfield(Par.NumericalParams, F),
        )
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

function setupDirectory(DirPath::String, Par; overwrite = false)
    DirPath = generateName_verbose(DirPath, Par)
    overwrite || (DirPath = UniqueDirName(DirPath))
    println("Checkpoints saved at $(abspath(DirPath))")
    # CheckPath = joinpath(DirPath,"Checkpoints")
    mkpath(DirPath)
    return DirPath
end

function setupDirectory(::Nothing, args...; kwargs...)
    return nothing
end

"""Rename CurrentState to FinalState as indicator that Job is finished"""
function SetCompletionCheckmark(DirPath::String)
    CurrState = joinpath(DirPath, "CurrentState.h5")
    FinState = UniqueFileName(joinpath(DirPath, "FinalState.h5"))
    if ispath(CurrState)
        mv(CurrState, FinState)
    end
end

SetCompletionCheckmark(::Nothing) = nothing

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
