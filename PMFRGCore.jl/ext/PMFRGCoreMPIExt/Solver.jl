using PMFRGCore, MPI, PencilArrays

function createLocalStateDecomposed(total_length::Int)::PencilArray

    decomposition = Pencil((total_length,), (1,), MPI.COMM_WORLD)

    PencilArray{Float64}(undef, decomposition)

end

# Utility to create a VBuffer out of a real vector knowing the pencil array geometry
# so that then MPI.AllgatherV! can be used to fill it in.
function getglobalbuff(
    globalState::Vector{T},
    localArray::PencilArray{T},
)::VBuffer where {T}
    nranks = length(localArray.pencil.topology.ranks)

    counts = [length(range_remote(localArray, rank)[1]) for rank = 1:nranks]


    VBuffer(globalState, counts)

end

# needs to create and ODE problem that uses the pencil arrays.
# the full arrays are passed inside Setup
function get_ODE_problem(globalState::Vector, t0, tend, setup, ps::PMFRGCore.UseMPI)
    Deriv_subst! = PMFRGCore.generateSubstituteDeriv(Deriv!, ps)
    State = createLocalStateDecomposed(length(globalState))
    ODEProblem(Deriv_subst!, State, (t0, tend), setup)
end



# Needs to create a derivative function that uses pencil arrays
function PMFRGCore.generateSubstituteDeriv(getDeriv!::Function, ::PMFRGCore.UseMPI)

    function DerivSubs!(Deriv::PencilArray, State::PencilArray, par, t)
        # println(t)
        Lam = t_to_Lam(t)
        (; globalDeriv, globalState) = par
        a = getDeriv!(globalDeriv, globalState, par, Lam)
        # Needs to copy the results into the pencil arrays
        Deriv .*= Lam
        a
    end

end

# Adds just the global state buffers (which are plain vectors)
# to the output of AllocateSetup(...,::PMFRGCore.MultiThreaded)
function PMFRGCore.AllocateSetup(Par::PMFRGCore.AbstractOneLoopParams, ::PMFRGCore.UseMPI)

    globalState = PMFRGCore.createStateVector(getArrayGeometry(Par); floattype = floattype)
    globalDeriv = similar(globalState)

    return merge(
        PMFRGCore.AllocateSetup(Par, PMFRGCore.MultiThreaded()),
        (; globalState, globalDeriv),
    )
end
