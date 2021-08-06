"""Saves Vertices to a compressed HDF5 file in a Group "Lam"."""
function saveState(Filename::String,State::ArrayPartition,Lam,mode = "cw")
    Vertices = State.x
    Names = "fint","gamma","Va","Vb","Vc"
    try
        h5open(Filename,mode) do f
            for (Name,V) in zip(Names,Vertices)
                f["States/$Lam/$Name",blosc = 9] = V
            end
            # h5write(Filename,"$Lam/$Name",compress(V,level = 9))
        end
        h5write(Filename,"States/$Lam/Lam",Lam)
    catch e
        @warn "Saving State was unsuccessfull with exception:\n $e "
    finally
        close(f)
    end
end

function getLowestKeyLambda(Filename::String)
    Lam_min = 1E10
    currkey = nothing
    h5open(Filename,"r") do f
        for k in keys(f["States"])
            Lam = Array(f["States/$k/Lam"])
            if Lam <Lam_min 
                Lam_min = Lam
                currkey = string("States/",k)
            end
        end
    end
    return currkey,Lam_min
end

"""Reads Vertices at lowest Lambda from file"""
function readLastState(Filename::String)
    Names = "fint","gamma","Va","Vb","Vc"
    lastkey,Lam = getLowestKeyLambda(Filename)
    State = ArrayPartition(
        (h5read(Filename,string(lastkey,"/$N")) for N in Names)...
    )
    return State,Lam
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
function saveParams(Filename,Par::Params)
    Fields = EssentialParamFields()
    for F in Fields
        h5write(Filename,"Params/$F",getfield(Par,F))
    end
    h5write(Filename,"Geometry/Name",Par.System.Name)
    h5write(Filename,"Geometry/NLen",Par.System.NLen)
    h5write(Filename,"Geometry/couplings",Par.System.couplings)
    h5write(Filename,"Geometry/Npairs",Par.System.Npairs)
end

function readParams(Filename,Geometry)
    Fields = EssentialParamFields()
    Kwargs = Dict((F => h5read(Filename,"Params/$F") for F in Fields)...)
    Par = Params(;System = Geometry,Kwargs...)
    return Par
end