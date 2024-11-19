using PMFRGCore, MPI, PencilArrays

function createGlobalStateVector(Par::PMFRGParams)

    PMFRGCore.createStateVector(
        PMFRGCore.getArrayGeometry(Par);
        floattype = PMFRGCore._getFloatType(Par),
    )
end

function PMFRGCore.InitializeState(Par::PMFRGParams, ::PMFRGCore.UseMPI)
    StateGlobalTmp = createGlobalStateVector(Par)

    Γc = PMFRGCore.getVc(StateGlobalTmp, PMFRGCore.getArrayGeometry(Par))

    (; couplings) = Par.System
    PMFRGCore.setToBareVertex!(Γc, couplings)

    decomposition = Pencil((length(StateGlobalTmp),), (1,), MPI.COMM_WORLD)

    floattype = PMFRGCore._getFloatType(Par)
    State = PencilArray{floattype}(undef, decomposition)
    State .= StateGlobalTmp[range_local(State)...]
    return State

end



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

    StateMPIBuff = VBuffer(
        let
            StateVector = createGlobalStateVector(Par)
            decomposition = Pencil((length(StateVector),), (1,), MPI.COMM_WORLD)
            counts = map(
                r -> length(range_remote(decomposition, r + 1)[1]),
                decomposition.topology.ranks,
            )
            StateVector, counts
        end...,
    )

    return merge(
        PMFRGCore.AllocateSetup(Par, PMFRGCore.MultiThreaded()),
        (; StateMPIBuff, ParallelizationScheme = UseMPI()),
    )
end
