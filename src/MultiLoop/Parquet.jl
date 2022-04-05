
"""Given set of parameters solve the self-consistent Parquet approximation iteratively. 
TODO: Save output to file"""
function SolveParquet(Par::ParquetParams,Lam::Real,getObsFunc::Function = getObservables;kwargs...)
    Workspace = SetupParquet(Par)
    SolveParquet(Workspace,Lam,getObsFunc;kwargs...)
end

function SolveParquet(Workspace::ParquetWorkspace,Lam::Real,getObsFunc::Function = getObservables;kwargs...)
    ObsType = typeof(getObservables(Workspace,Lam))
    Obs = StructArray(ObsType[])
    @time iterateSolution!(Workspace,Lam,Obs,getObsFunc)
end

"""Obtains a solution to Bethe-Salpeter and Schwinger-Dyson equations by iteration until convergence is reached up to accuracy specified by accuracy in Params"""
function iterateSolution!(Workspace::ParquetWorkspace,Lam::Real,Obs,getObsFunc::Function)
    @unpack OldState,State,I,Γ0,X,B0,BX,Par,Buffer = Workspace
    
    maxIterBSE = Par.Options.maxIterBSE


    Tol_Vertex = 1E16 #Initial tolerance level

    iG(x,nw) = iG_(State.γ,x,Lam,nw,T)

    getProp! = constructPropagatorFunction(Workspace,Lam)
    iter = 0
    while Tol_Vertex > Par.NumericalParams.accuracy 
        iter+=1
        
        if iter > maxIterBSE
            @warn("BSE: No convergence found after $maxIterBSE iterations")
            break
        end

        writeTo!(OldState.Γ,State.Γ)
        
        # iterateSDE!(Workspace,Lam)
        
        computeLeft2PartBubble!(B0,Γ0,Γ0,State.Γ,getProp!,Par,Buffer)
        computeLeft2PartBubble!(BX,X,X,State.Γ,getProp!,Par,Buffer)
        
        getXFromBubbles!(X,B0,BX) #TODO when applicable, this needs to be generalized for beyond-Parquet approximations 
        getVertexFromChannels!(State.Γ,I,X)
        symmetrizeVertex!(State.Γ,Par)
        
        iterateSDE!(Workspace,Lam)
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

function constructPropagatorFunction(Workspace::PMFRGWorkspace,Lam)    
    Par = Workspace.Par
    T= Par.NumericalParams.T
    NUnique = Par.System.NUnique

    @inline iG(x,nw) = iG_(Workspace.State.γ,x,Lam,nw,T)

    function getProp!(BubbleProp,nw1,nw2)
        for i in 1:NUnique, j in 1:NUnique
            BubbleProp[i,j] = iG(i,nw1) *iG(j,nw2)* Par.NumericalParams.T
        end
        return BubbleProp
    end
    return getProp!
end

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


function iterateSDE!(Workspace::ParquetWorkspace,Lam)
    @unpack OldState,State,Γ0,X,B0,BX,Par,Buffer = Workspace
    @inline Prop(x,nw) = 1/6*iG_(OldState.γ,x,Lam,nw,Par.NumericalParams.T)#*3

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

        SDE_tolerance = reldist(State.γ,OldState.γ)
    end
    if !Par.Options.MinimalOutput
        println("""
        \t\tSDE step done after $iter / $maxIterSDE iterations (tol = $(SDE_tolerance)""")
    end
    return 
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
