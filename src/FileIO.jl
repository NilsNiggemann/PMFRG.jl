"""
This file contains functions used for reading and writing FRG data to disc.
Stored data can be used to re-launch an incomplete FRG calculation.
"""

"""Saves Vertices to a compressed HDF5 file in a Group "Lam"."""
function saveState(Filename::String,State::ArrayPartition,Lam,mode = "cw")
    Vertices = State.x
    Names = "fint","gamma","Va","Vb","Vc"
    # Filename = string(DirName,"/$(string(round(Lam,digits =3))).h5")
    try
        h5open(Filename,mode) do f
            for (Name,V) in zip(Names,Vertices)
                f["$Name",blosc = 9] = V
            end
        end
        h5write(Filename,"Lam",Lam)
    catch e
        @warn "Saving state was unsuccessfull with exception:\n $e "
    end
end

"""Reads Vertices from file"""
function readState(Filename::String)
    Names = "fint","gamma","Va","Vb","Vc"
    State = ArrayPartition(
        (h5read(Filename,"$N") for N in Names)...
    )
    return State
end

function readLam(Filename::String)
    Lam = h5read(Filename,"Lam")
    return Lam
end

function readObservables(Filename::String)
    t = h5read(Filename,"Observables/Lambda")
    Fields = fieldnames(Observables)
    obsTuple = Tuple(h5read(Filename,"Observables/$f") for f in Fields)
    saveval = Observables[]
    saved_values = SavedValues(double,Observables)
    # return obsTuple
    for i in eachindex(t)
        Currentval(x) = selectdim(x,length(size(x)),i) |>Array
        # return Currentval.(obsTuple)
        push!(saveval,Observables(Currentval.(obsTuple)...))
    end
    append!(saved_values.t,t)
    append!(saved_values.saveval, saveval)
    return saved_values
end

function readGeometry(Filename::String, GeometryGenerator::Function;Group="Geometry")
    NLen = h5read(Filename,"$Group/NLen")
    couplings = h5read(Filename,"$Group/couplings")
    Npairs = h5read(Filename,"$Group/Npairs")
    Name = h5read(Filename,"$Group/Name")
    System = GeometryGenerator(NLen)
    @assert Npairs == System.Npairs "Number of unique pairs does not match geometry in $Filename"
    @assert Name == System.Name "Name does not match geometry in $Filename"
    System.couplings .= couplings
    return System
end

EssentialParamFields() = (
    :T,
    :N,
    :Ngamma,
    :accuracy,
    :Lam_min,
    :Lam_max,
    :ex_freq,
    :lenIntw,
    :lenIntw_acc
)

"""Saves important information about computation parameters so that they can be reconstructed"""
function saveParams(Filename,Par::Params,Group = "")
    Fields = EssentialParamFields()
    for F in Fields
        h5write(Filename,joinpath(Group,"Params/$F"),getfield(Par,F))
    end
    h5write(Filename,joinpath(Group,"Geometry/Name"),Par.System.Name)
    h5write(Filename,joinpath(Group,"Geometry/NLen"),Par.System.NLen)
    h5write(Filename,joinpath(Group,"Geometry/couplings"),Par.System.couplings)
    h5write(Filename,joinpath(Group,"Geometry/Npairs"),Par.System.Npairs)
end

function readParams(Filename::String,Geometry::SpinFRGLattices.Geometry;modifyParams...)
    Fields = EssentialParamFields()
    Kwargs = Dict((F => h5read(Filename,"Params/$F") for F in Fields)...)
    Par = Params(;System = Geometry,Kwargs...,modifyParams...)
    return Par
end

readParams(Filename::String,GeometryGenerator::Function;modifyParams...) = readParams(Filename,readGeometry(Filename,GeometryGenerator);modifyParams...)

function modifyParams(Par;modifyParams...)
    ParKwargs = Dict((F => getfield(Par,F) for F in fieldnames(Params))...)
    Par = Params(;ParKwargs...,modifyParams...)
end

function setupDirectory(DirPath,Par)
    DirPath = generateUniqueName(DirPath,Par)
    println("Checkpoints saved at $DirPath")
    # CheckPath = joinpath(DirPath,"Checkpoints")
    mkpath(DirPath)
    return DirPath
end

function saveCurrentState(DirPath::String,State::AbstractArray,saved_Values::DiffEqCallbacks.SavedValues,Lam::Real,Par::Params)
    saveState(joinpath(DirPath,"CurrentState.h5"),State,Lam,"w")
    saveParams(joinpath(DirPath,"CurrentState.h5"),Par)
    saveObs(joinpath(DirPath,"CurrentState.h5"),saved_Values,"Observables")
end

"""Rename CurrentState to FinalState as indicator that Job is finished"""
function SetCompletionCheckmark(DirPath::String)
    CurrState = joinpath(DirPath,"CurrentState.h5")
    FinState = joinpath(DirPath,"FinalState.h5")
    if ispath(CurrState) && !ispath(FinState)
        mv(joinpath(DirPath,"CurrentState.h5"),joinpath(DirPath,"FinalState.h5"))
    end
end

function UniqueDirName(FullPath)
    newpath = FullPath
    versionPath(index) = string(FullPath,"(v_",index,")")
    while ispath(newpath)
        it = findfirst("(v_",newpath)
        if it !== nothing
            currindex = parse(Int,newpath[it[end]+1:end-1])
            newpath =versionPath(currindex+1)
        else
            newpath =versionPath(1)
        end
    end
    return newpath
end


function generateUniqueName(Directory::String,Par::Params)
    @unpack Name,T,N = Par
    Name = "$(Name)_N=$(N)_T=$T"
    Name = UniqueDirName(joinpath(Directory,Name))
    return Name
end

function generateFileName(Par::Params,arg::String = "")
    @unpack Name,N = Par
    Name = "$(Name)_N=$(N)_l_1_$arg.h5"
    return Name
end

function ParamsCompatible(Par1,Par2)
    Fields = (:Npairs,:N,:Ngamma,:Lam_max,:System) 
    for f in Fields
        @assert getproperty(Par1,f) == getproperty(Par2,f) "$f not compatible with parameters used in Checkpoint"
    end

    return true
end

function getFileParams(Filename,Geometry,Par::Params)
    Lam = readLam(Filename)
    Par_file = readParams(Filename,Geometry;Lam_max = Lam)
    Par = modifyParams(Par;Lam_max = Lam)
    ParamsCompatible(Par,Par_file)
    return Par
end

function getFileParams(Filename,Geometry,Par::Nothing)
    Lam = readLam(Filename)
    return readParams(Filename,Geometry;Lam_max = Lam)
end

function SolveFRG_Checkpoint(Filename::String,Geometry::SpinFRGLattices.Geometry,Par = nothing;MainFile = nothing,Group =nothing,overwrite=false,kwargs...)
    State = readState(Filename)
    Old_Lam_max = h5read(Filename,"Params/Lam_max")
    Par = getFileParams(Filename,Geometry,Par)
    saved_values_full = readObservables(Filename)
    setup = AllocateSetup(Par)
    CheckPointfolder = dirname(Filename)
    overwrite && mv(CheckPointfolder,CheckPointfolder*"_OLD")

    FilePath = dirname(CheckPointfolder)
    ObsSaveat = getLambdaMesh(nothing,Par.Lam_min,Old_Lam_max)
    filter!(x-> x<Par.Lam_max,ObsSaveat)
    sol,saved_values = launchPMFRG!(State,setup,getDeriv!;CheckpointDirectory = FilePath,ObsSaveat = ObsSaveat, kwargs...)
    append!(saved_values_full.t,saved_values.t)
    append!(saved_values_full.saveval,saved_values.saveval)
    if MainFile !== nothing
        saveMainOutput(MainFile,sol,saved_values_full,Par,Group)
    end
    overwrite && rm(CheckPointfolder*"_OLD",recursive=true)
    return sol,saved_values_full
end

function SolveFRG_Checkpoint(Filename::String,GeometryGenerator::Function,Par = nothing;kwargs...)
    System = readGeometry(Filename,GeometryGenerator)
    SolveFRG_Checkpoint(Filename,System,Par;kwargs...)
end

"""Saves Observables"""
function saveObs(Filename,saved_values::DiffEqCallbacks.SavedValues,Group = "")
    Fields = fieldnames(eltype(saved_values.saveval))
    ObsArr = StructArray(saved_values.saveval)
    for F in Fields
        arr = convertToArray(getproperty(ObsArr,F))
        h5write(Filename,joinpath(Group,string(F)),arr)
    end
    h5write(Filename,joinpath(Group,"Lambda"),saved_values.t)
end

function convertToArray(VecOfArray::AbstractVector{VT}) where {T,N,VT <: AbstractArray{T,N}}
    cat(VecOfArray...,dims = N+1)
end

function saveMainOutput(Filename::String,saved_values::DiffEqCallbacks.SavedValues,Group::String)
    println("Saving Main output to ", Filename)
    mkpath(dirname(Filename))
    saveObs(Filename,saved_values,Group)
end

function saveMainOutput(Filename::String,Solution::ODESolution,saved_values::DiffEqCallbacks.SavedValues,Par::Params,Group::String)
    @unpack T,System,N = Par
    @unpack Name,NUnique,Npairs,NLen = System
    Lambda = saved_values.t
    function save(name)
        saveMainOutput(name,saved_values,Group)
        Chi_nu = getChi(Solution[end],Lambda[end],Par,N)
        h5write(name,"$Group/T",T)
        h5write(name,"$Group/N",N)
        h5write(name,"$Group/NUnique",NUnique)
        h5write(name,"$Group/NLen",NLen)
        h5write(name,"$Group/Chi_nu",Chi_nu)
    end
    try
        save(Filename)
    catch e
        newName = UniqueDirName(Filename)
        @warn "Writing to $Filename errored with exception $(string(e))! Writing to $newName instead."
        save(newName)
    end
    # saveParams(Filename,Par,Group)
end

saveMainOutput(Filename::String,Solution::ODESolution,saved_values::DiffEqCallbacks.SavedValues,Par::Params,Group::Nothing) = saveMainOutput(Filename,Solution,saved_values,Par,string(Par.T))

function getFilesFromSubDirs(Folder::String)
    allpaths = collect(walkdir(Folder))
    Files = String[]
    for p in allpaths
        for filename in p[end]
            pathAndName = joinpath(p[begin],filename)
            push!(Files,pathAndName)
        end
    end
    return Files
end

function getUnfinishedJobs(Folder::String)
    allFiles = filter!(x->occursin("CurrentState.h5",x),getFilesFromSubDirs(Folder))

    getLam_min(file) = h5read(file,"Params/Lam_min")
    filter!(x-> getLam_min(x) != readLam(x),allFiles)
end