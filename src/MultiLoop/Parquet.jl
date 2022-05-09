
"""Given set of parameters solve the self-consistent Parquet approximation iteratively. 
TODO: Save output to file"""
function SolveParquet(Par::ParquetParams,Lam::Real,getObsFunc::Function = getObservables;kwargs...)
    Workspace = SetupParquet(Par)
    SolveParquet(Workspace,Lam,getObsFunc;kwargs...)
end

function SolveParquet(Workspace::ParquetWorkspace,Lam::Real,getObsFunc::Function = getObservables;kwargs...)
    ObsType = typeof(getObservables(Workspace,Lam))
    Obs = StructArray(ObsType[])
    @time iterateSolution_FP!(Workspace,Lam,Obs,getObsFunc)
    # @time iterateSolution!(Workspace,Lam,Obs,getObsFunc)
end

"""Obtains a solution to Bethe-Salpeter and Schwinger-Dyson equations by iteration until convergence is reached up to accuracy specified by accuracy in Params"""
function iterateSolution!(Workspace::ParquetWorkspace,Lam::Real,Obs,getObsFunc::Function;SDECheatfactor = 1)
    @unpack OldState,State,I,Γ0,X,B0,BX,Par,Buffer = Workspace
    
    maxIterBSE = Par.Options.maxIterBSE


    Tol_Vertex = 1E16 #Initial tolerance level

    getProp! = constructPropagatorFunction(Workspace,Lam)
    iter = 0
    while Tol_Vertex > Par.NumericalParams.accuracy 
        iter+=1
        
        if iter > maxIterBSE
            @warn("BSE: No convergence found after $maxIterBSE iterations, Tol = $Tol_Vertex")
            break
        end

        writeTo!(OldState.Γ,State.Γ)
        
        # iterateSDE!(Workspace,Lam)
        
        computeLeft2PartBubble!(B0,Γ0,Γ0,State.Γ,getProp!,Par,Buffer)
        computeLeft2PartBubble!(BX,X,X,State.Γ,getProp!,Par,Buffer)
        
        getXFromBubbles!(X,B0,BX) #TODO when applicable, this needs to be generalized for beyond-Parquet approximations 
        getVertexFromChannels!(State.Γ,I,X)
        symmetrizeVertex!(State.Γ,Par)
        
        # iterateSDE!(Workspace,Lam,SDECheatfactor = SDECheatfactor)
        iterateSDE_FP!(Workspace,Lam)
        Tol_Vertex = reldist(OldState.Γ,State.Γ)
        
        CurrentObs = getObsFunc(Workspace,Lam)
        push!(Obs,CurrentObs)

        if !Par.Options.MinimalOutput
            println("iteration $iter:")
            writeOutput(State,CurrentObs,Lam,Par)
            println("""
            Tol_Vertex = $Tol_Vertex
            """)
        end
        isnan(Tol_Vertex) && @warn ": BSE: Solution diverged after $iter iterations\n"
    end
    return Workspace,Obs
end


function constructPropagatorFunction(γ, Lam,Par)    
    T= Par.NumericalParams.T
    NUnique = Par.System.NUnique

    @inline iG(x,nw) = iG_(γ,x,Lam,nw,T)

    function getProp!(BubbleProp,nw1,nw2)
        for i in 1:NUnique, j in 1:NUnique
            BubbleProp[i,j] = iG(i,nw1) *iG(j,nw2)* T
        end
        return BubbleProp
    end
    return getProp!
end

constructPropagatorFunction(Workspace::PMFRGWorkspace,Lam) = constructPropagatorFunction(Workspace.State.γ,Lam,Workspace.Par) 




function getXFromBubbles!(X::BubbleType,B0::BubbleType,BX::BubbleType)
    for f in fieldnames(BubbleType)
        Xf,B0f,BXf = getfield(X,f),getfield(B0,f),getfield(BX,f)

        Threads.@threads for i in eachindex(Xf,B0f,BXf)
           Xf[i] = 0.5* B0f[i] + BXf[i]
        end
    end
    return X
end

function getVertexFromChannels!(Γ::VertexType,I::VertexType,X::BubbleType)
    Threads.@threads for iu in axes(Γ.a,4)
        for it in axes(Γ.a,3), is in axes(Γ.a,2), Rij in axes(Γ.a,1)
            Γ.a[Rij,is,it,iu] = I.a[Rij,is,it,iu] + X.a[Rij,is,it,iu] - X.Ta[Rij,it,is,iu] + X.Ta[Rij,iu,is,it]
            Γ.b[Rij,is,it,iu] = I.b[Rij,is,it,iu] + X.b[Rij,is,it,iu] - X.Tc[Rij,it,is,iu] + X.Tc[Rij,iu,is,it]
            Γ.c[Rij,is,it,iu] = I.c[Rij,is,it,iu] + X.c[Rij,is,it,iu] - X.Tb[Rij,it,is,iu] + X.Td[Rij,iu,is,it]
        end
    end
    return Γ
end

"""Self-consistently iterates SDE until convergence is reached. SDECheatfactor only to be used for debugging purposes. Setting it to 3 seems to reproduce PMFRG for some reason."""
function iterateSDE!(Workspace::ParquetWorkspace,Lam;SDECheatfactor = 1)
    @unpack OldState,State,Γ0,X,B0,BX,Par,Buffer = Workspace
    @inline Prop(x,nw) = 1/6*iG_(OldState.γ,x,Lam,nw,Par.NumericalParams.T)*SDECheatfactor

    getProp! = constructPropagatorFunction(Workspace,Lam)

    maxIterSDE = Par.Options.maxIterSDE
    SDE_tolerance = 1E16
    iter = 0
    while SDE_tolerance > Par.Options.SDE_tolerance
        if iter >= maxIterSDE
            # @warn("SDE: No convergence found after $maxIter iterations\nremaining Tol: $SDE_tolerance")
            break
        end
        iter +=1
        writeTo!(OldState.γ,State.γ)

        # computeLeft2PartBubble!(B0,Γ0,Γ0,State.Γ,getProp!,Par,Buffer)
        # println(getProp!(Buffer.Props[1],1,1))
        # iter < 10 && println(State.γ[1,1:5])
        
        compute1PartBubble!(State.γ,B0,Prop,Par)
        # compute1PartBubble_BS!(State.γ,State.Γ,Γ0,Prop,Par)

        SDE_tolerance = reldist(State.γ,OldState.γ)
    end
    if !Par.Options.MinimalOutput
        println("""
        \t\tSDE step done after $iter / $maxIterSDE iterations (tol = $(SDE_tolerance)""")
    end
    return 
end
function iterateSolution_FP!(Workspace::ParquetWorkspace,Lam::Real,Obs,getObsFunc::Function)
    @unpack OldState,State,I,Γ0,X,B0,BX,Par,Buffer = Workspace
    T = Par.NumericalParams.T
    maxIterBSE = Par.Options.maxIterBSE

    OldStateArr,StateArr = ArrayPartition.((OldState,State))
        
    function FixedPointFunction!(State_Arr,OldState_Arr)
        anyisnan(OldState_Arr) && return State_Arr
        gamma = OldState_Arr.x[2]

        State = StateType(State_Arr.x...)
        OldState = StateType(OldState_Arr.x...)
        getProp! = constructPropagatorFunction(gamma,Lam,Par)

        computeLeft2PartBubble!(B0,Γ0,Γ0,OldState.Γ,getProp!,Par,Buffer)
        computeLeft2PartBubble!(BX,X,X,OldState.Γ,getProp!,Par,Buffer)
        
        getXFromBubbles!(X,B0,BX) #TODO when applicable, this needs to be generalized for beyond-Parquet approximations 
        getVertexFromChannels!(State.Γ,I,X)
        symmetrizeVertex!(State.Γ,Par)
        
        
        # iterateSDE_FP!(Workspace,Lam) #Todo: underlying function takes Old state as an input, but not the new state computed here?
        # iterateSDE_FP!(State.γ,State.Γ,Γ0,Lam,Par)
        iterateSDE_FP!(State.γ,State.Γ,B0,Γ0,Lam,Par,Buffer)
        WS_State = ArrayPartition(Workspace.State)
        WS_State .= State_Arr
        CurrentObs = getObsFunc(Workspace,Lam)
        push!(Obs,CurrentObs)
        if !Par.Options.MinimalOutput
            writeOutput(State,CurrentObs,Lam,Par)
        end
        # OldState_Arr .= State_Arr
        return State_Arr
    end
    
    s = afps!(FixedPointFunction!,OldStateArr,iters = maxIterBSE,vel = 0.0,ep = getEpsilon(T),tol = Par.Options.SDE_tolerance)
    StateArr .= s.x
    println("""
    \t\tBSE done after  $(s.iters) / $maxIterBSE iterations (tol = $(s.error))""")
    anyisnan(StateArr) && @warn "NaN detected... Aborted"
    return Workspace,Obs
end

anyisnan(A::ArrayPartition) = any((any(isnan,i) for i in A.x))


function iterateSDE_FP!(γ,Γ,B0,Γ0,Lam,Par,Buffer)
    maxIterSDE = Par.Options.maxIterSDE
    T = Par.NumericalParams.T
    function FixedPointFunction!(gamma,gammaOld)
        @inline Prop(x,nw) = 1/6*iG_(gammaOld,x,Lam,nw,T)#*3

        getProp! = constructPropagatorFunction(gammaOld,Lam,Par)

        # computeLeft2PartBubble!(B0,Γ0,Γ0,Γ,getProp!,Par,Buffer)
        compute1PartBubble!(gamma,B0,Prop,Par)
        # compute1PartBubble_BS!(gamma,Γ,Γ0,Prop,Par)
        return gamma
    end

    s = afps!(FixedPointFunction!,γ,iters = maxIterSDE,vel = 0.2,ep = getEpsilon(T),tol = Par.Options.SDE_tolerance)
    γ .= s.x
    # FixedPointFunction!(State.γ,OldState.γ)
    # println(maximum(γ))
    # println(maximum(ArrayPartition(Γ.a,Γ.b,Γ.c)))
    if !Par.Options.MinimalOutput
        println("""
        \t\tSDE step done after  $(s.iters) / $maxIterSDE iterations (tol = $(s.error))""")
    end
end

function getEpsilon(T)
    min(T/2,1.)
end

getChi(State::StateType, Lam::Real,Par::ParquetParams) = getChi(State.γ,State.Γ.c, Lam,Par)
getChi(State::StateType, Lam::Real,Par::ParquetParams,Numax) = getChi(State.γ,State.Γ.c, Lam,Par,Numax)

function getObservables(Workspace::ParquetWorkspace,Lam)
    State = Workspace.State
    chi = getChi(Workspace.State,Lam,Workspace.Par)
    MaxVa = maximum(abs,State.Γ.a,dims = (2,3,4,5))[:,1,1,1]
    MaxVb = maximum(abs,State.Γ.b,dims = (2,3,4,5))[:,1,1,1]
    MaxVc = maximum(abs,State.Γ.c,dims = (2,3,4,5))[:,1,1,1]
    return Observables(chi,copy(State.γ),copy(State.f_int),MaxVa,MaxVb,MaxVc) # make sure to allocate new memory each time this function is called
end

# writeOutput(St::StateType,Obs,Lam,Par) = println(Obs)
writeOutput(St::StateType,Obs,Lam,Par) = writeOutput(St.f_int,St.γ,St.Γ.a,St.Γ.b,St.Γ.c,Obs,Lam,Par)

function modifyWorkspace(Workspace::ParquetWorkspace;kwargs...)
    NewPar = modifyParams(Workspace.Par;kwargs...)
    newWorkspace = SetupParquet(NewPar)
    writeTo!(newWorkspace.State,Workspace.State)
    return newWorkspace
end
