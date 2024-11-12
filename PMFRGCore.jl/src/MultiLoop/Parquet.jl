
"""Performs single BSE iteration. Specify State explicitly to make fixed point libraries compatible (even though it may point to the same memory location as Workspace.State) """
function BSE_iteration!(State::StateType, Workspace::ParquetWorkspace, Lam::Real)
    (; OldState, I, Γ0, X, B0, BX, Par, Buffer) = Workspace

    getProp! = constructPropagatorFunction(OldState.γ, Lam, Par)
    computeLeft2PartBubble!(B0, Γ0, Γ0, OldState.Γ, getProp!, Par, Buffer)
    computeLeft2PartBubble!(BX, X, X, OldState.Γ, getProp!, Par, Buffer)

    getXFromBubbles!(X, B0, BX) #TODO when applicable, this needs to be generalized for beyond-Parquet approximations 
    getVertexFromChannels!(State.Γ, I, X)
    symmetrizeVertex!(State.Γ, Par)


    return State
end

"""Obtains a solution to Bethe-Salpeter and Schwinger-Dyson equations by iteration until convergence is reached up to accuracy specified by accuracy in Params"""
function iterateSolution!(Workspace::ParquetWorkspace, Lam::Real, Obs, getObsFunc::Function)
    (; OldState, State, Par) = Workspace

    BSE_iters = Par.Options.BSE_iters


    Tol_Vertex = 1E16 #Initial tolerance level

    iter = 0
    while Tol_Vertex > Par.NumericalParams.accuracy
        iter += 1

        if iter >= BSE_iters
            @warn(
                "BSE: No convergence found after $BSE_iters iterations, tol = $Tol_Vertex"
            )
            break
        end

        writeTo!(OldState.Γ, State.Γ)
        BSE_iteration!(Workspace.State, Workspace, Lam)
        # iterateSDE_FP!(State.γ,State.Γ,B0,Γ0,Lam,Par,Buffer)
        Tol_Vertex = reldist(OldState.Γ, State.Γ)

        iterateSDE!(Workspace, Lam)
        dampen!(State, OldState, Par.Options.BSE_epsilon)

        CurrentObs = getObsFunc(Workspace, Lam)
        push!(Obs, CurrentObs)

        if !Par.Options.MinimalOutput
            println("iteration $iter:")
            writeOutput(State, CurrentObs, Lam, Par)
            println("""
            Tol_Vertex = $Tol_Vertex
            """)
        end
        isnan(Tol_Vertex) && @warn ": BSE: Solution diverged after $iter iterations\n"
    end
    Par.Options.MinimalOutput || println("""
    \t\tBSE done after  $(iter) / $BSE_iters iterations (tol = $(Tol_Vertex))""")

    return Workspace, Obs
end

function dampen!(StateArr::AbstractArray, OldStateArr::AbstractArray, epsilon::Real)
    @. StateArr = (1 - epsilon) * OldStateArr + epsilon * StateArr
end

function constructPropagatorFunction(γ, Lam, Par)
    T = Par.NumericalParams.T
    NUnique = Par.System.NUnique

    @inline iG(x, nw) = iG_(γ, x, Lam, nw, T)

    function getProp!(BubbleProp, nw1, nw2)
        for i = 1:NUnique, j = 1:NUnique
            BubbleProp[i, j] = iG(i, nw1) * iG(j, nw2) * T
        end
        return BubbleProp
    end
    return getProp!
end

constructPropagatorFunction(Workspace::PMFRGWorkspace, Lam) =
    constructPropagatorFunction(Workspace.State.γ, Lam, Workspace.Par)


function getXFromBubbles!(X::BubbleType, B0::BubbleType, BX::BubbleType)
    for f in fieldnames(BubbleType)
        Xf, B0f, BXf = getfield(X, f), getfield(B0, f), getfield(BX, f)

        Threads.@threads for i in eachindex(Xf, B0f, BXf)
            Xf[i] = 0.5 * B0f[i] + BXf[i]
        end
    end
    return X
end

function getVertexFromChannels!(Γ::VertexType, I::VertexType, X::BubbleType)
    Threads.@threads for iu in axes(Γ.a, 4)
        for it in axes(Γ.a, 3), is in axes(Γ.a, 2), Rij in axes(Γ.a, 1)
            Γ.a[Rij, is, it, iu] =
                I.a[Rij, is, it, iu] + X.a[Rij, is, it, iu] - X.Ta[Rij, it, is, iu] +
                X.Ta[Rij, iu, is, it]
            Γ.b[Rij, is, it, iu] =
                I.b[Rij, is, it, iu] + X.b[Rij, is, it, iu] - X.Tc[Rij, it, is, iu] +
                X.Tc[Rij, iu, is, it]
            Γ.c[Rij, is, it, iu] =
                I.c[Rij, is, it, iu] + X.c[Rij, is, it, iu] - X.Tb[Rij, it, is, iu] +
                X.Td[Rij, iu, is, it]
        end
    end
    return Γ
end

"""Self-consistently iterates SDE until convergence is reached."""
function iterateSDE!(Workspace::ParquetWorkspace, Lam)
    (; OldState, State, Γ0, X, B0, BX, Par, Buffer) = Workspace
    @inline Prop(x, nw) = 1 / 6 * iG_(OldState.γ, x, Lam, nw, Par.NumericalParams.T)

    # getProp! = constructPropagatorFunction(Workspace,Lam)

    SDE_iters = Par.Options.SDE_iters
    SDE_tolerance = 1E16
    iter = 0
    while SDE_tolerance > Par.Options.SDE_tolerance
        if iter >= SDE_iters
            # @warn("SDE: No convergence found after $maxIter iterations\nremaining Tol: $SDE_tolerance")
            break
        end
        iter += 1
        writeTo!(OldState.γ, State.γ)

        # computeLeft2PartBubble!(B0,Γ0,Γ0,State.Γ,getProp!,Par,Buffer)

        compute1PartBubble!(State.γ, B0, Prop, Par)
        # compute1PartBubble_BS!(State.γ,State.Γ,Γ0,Prop,Par)
        dampen!(State.γ, OldState.γ, Par.Options.SDE_epsilon)
        SDE_tolerance = reldist(State.γ, OldState.γ)
    end
    if !Par.Options.MinimalOutput
        println("""
        \t\tSDE step done after $iter / $SDE_iters iterations (tol = $(SDE_tolerance)""")
    end
    return
end



function iterateSolution_FP!(Workspace::ParquetWorkspace, Lam::Real, Obs)
    (; OldState, State, Γ0, B0, Par, Buffer) = Workspace
    (; BSE_iters, BSE_epsilon, BSE_vel) = Par.Options
    (; accuracy) = Par.NumericalParams

    ObsType = eltype(Obs)
    OldStateArr, StateArr = repackStateVector.((OldState, State))

    function FixedPointFunction!(State_Arr, OldState_Arr)
        any(isnan, OldState_Arr) && return State_Arr
        unpackStateVector!(State, State_Arr)
        BSE_iteration!(State, Workspace, Lam)
        iterateSDE_FP!(State.γ, State.Γ, B0, Γ0, Lam, Par, Buffer)

        writeTo!(Workspace.State, State) # need to do this since State is not equal to Workspace.State anymore!

        CurrentObs = getObservables(ObsType, Workspace, Lam)
        push!(Obs, CurrentObs)
        if !Par.Options.MinimalOutput
            writeOutput(State, CurrentObs, Lam, Par)
        end
        # OldState_Arr .= State_Arr
        repackStateVector!(State_Arr, State)
        return State_Arr
    end

    s = afps!(
        FixedPointFunction!,
        OldStateArr,
        iters = BSE_iters,
        vel = BSE_vel,
        ep = BSE_epsilon,
        tol = accuracy,
    )
    StateArr .= s.x
    println("""
    \t\tBSE done after  $(length(Obs)) / $BSE_iters iterations (tol = $(s.error))""")
    if any(isnan, StateArr)
        @warn "NaN detected... Aborted"
    elseif s.error > accuracy
        @warn "Tolerance goal ($accuracy) not reached"
    end
    return Workspace, Obs
end



function iterateSDE_FP!(γ, Γ, B0, Γ0, Lam, Par, Buffer)
    (; SDE_iters, SDE_vel, SDE_epsilon) = Par.Options
    T = Par.NumericalParams.T
    function FixedPointFunction!(gamma, gammaOld)
        @inline Prop(x, nw) = 1 / 6 * iG_(gammaOld, x, Lam, nw, T)

        getProp! = constructPropagatorFunction(gammaOld, Lam, Par)

        # computeLeft2PartBubble!(B0,Γ0,Γ0,Γ,getProp!,Par,Buffer)
        compute1PartBubble!(gamma, B0, Prop, Par)
        # compute1PartBubble_BS!(gamma,Γ,Γ0,Prop,Par)
        return gamma
    end

    s = afps!(
        FixedPointFunction!,
        γ,
        iters = SDE_iters,
        vel = SDE_vel,
        ep = SDE_epsilon,
        tol = Par.Options.SDE_tolerance,
    )
    γ .= s.x
    # FixedPointFunction!(State.γ,OldState.γ)
    # println(maximum(γ))
    if !Par.Options.MinimalOutput
        println("""
        \t\tSDE step done after  $(s.iters) / $SDE_iters iterations (tol = $(s.error))""")
    end
end


getChi(State::StateType, Lam::Real, Par::ParquetParams) =
    getChi(State.γ, State.Γ.c, Lam, Par)
getChi(State::StateType, Lam::Real, Par::ParquetParams, Numax) =
    getChi(State.γ, State.Γ.c, Lam, Par, Numax)

function getObservables(::Type{Observables}, Workspace::ParquetWorkspace, Lam)
    State = Workspace.State
    chi = getChi(Workspace.State, Lam, Workspace.Par)
    MaxVa = maximum(abs, State.Γ.a, dims = (2, 3, 4, 5))[:, 1, 1, 1]
    MaxVb = maximum(abs, State.Γ.b, dims = (2, 3, 4, 5))[:, 1, 1, 1]
    MaxVc = maximum(abs, State.Γ.c, dims = (2, 3, 4, 5))[:, 1, 1, 1]
    return Observables(chi, copy(State.γ), copy(State.f_int), MaxVa, MaxVb, MaxVc) # make sure to allocate new memory each time this function is called
end

# writeOutput(St::StateType,Obs,Lam,Par) = println(Obs)
writeOutput(St::StateType, Obs, Lam, Par) =
    writeOutput(St.f_int, St.γ, St.Γ.a, St.Γ.b, St.Γ.c, Obs, Lam, Par)

function modifyWorkspace(Workspace::ParquetWorkspace; kwargs...)
    NewPar = modifyParams(Workspace.Par; kwargs...)
    newWorkspace = SetupParquet(NewPar)
    writeTo!(newWorkspace.State, Workspace.State)
    return newWorkspace
end
