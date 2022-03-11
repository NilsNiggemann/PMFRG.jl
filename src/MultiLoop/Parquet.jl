struct ParquetWorkspace{T,P <: Params{S,M} where {S,M}}
    State::StateType{T} #Stores the current step of iteration
    OldState::StateType{T} #Stores the last step of iteration

    Γ0::BareVertexType{T} #Stores the bare vertex
    I::VertexType{T} #Stores the irreducible vertex
    
    X::BubbleType{T} #Stores the actual bubble funtion X = B0 + BX

    B0::BubbleType{T} # Bubble involving bare vertex Γ_0 ∘ Γ
    BX::BubbleType{T} # Bubble involving other bubble X ∘ Γ

    Par::P # Params
end

"""Constructs Workspace for parquet equations. Irreducible vertex can later be provided by optional argument but this is currently not supported """
function SetupParquet(Par)
    IrrVertex = BareVertex_Freq(Par) 
    Workspace = ParquetWorkspace(
        StateType(Par),
        StateType(Par),

        BareVertexType(Par),
        IrrVertex,
        
        BubbleType(Par),

        BubbleType(Par),
        BubbleType(Par),
        
        Par
        )
    setToBareVertex!(Workspace.OldState.Γ,Par)
    return Workspace
end

"""Given set of parameters solve the self-consistent Parquet approximation iteratively. 
TODO: Save output to file"""
function SolveParquet(Par::Params,Lam::Real,getObsFunc::Function = getObservables;maxiterBubble::Integer = 10, maxitergamma::Integer=100,kwargs...)
    Workspace = SetupParquet(Par)
    
    PropsBuffers = [Matrix{double}(undef,Par.NUnique,Par.NUnique) for _ in 1:Threads.nthreads()] 
    VertexBuffers = [VertexBuffer(Par.Npairs) for _ in 1:Threads.nthreads()] 
    
    ObsType = typeof(getObsFunc(Workspace,Lam))
    Obs = StructArray(ObsType[])

    @time iterateSolution!(Workspace,Lam,PropsBuffers,VertexBuffers,maxiterBubble,maxitergamma,Obs,getObsFunc)
end

"""Obtains a solution to Bethe-Salpeter and Schwinger-Dyson equations by iteration until convergence is reached up to accuracy specified by accuracy in Params"""
function iterateSolution!(Workspace::ParquetWorkspace,Lam::Real,PropsBuffers,VertexBuffers,maxiterBubble::Integer, maxitergamma::Integer,Obs,getObsFunc::Function)
    @unpack OldState,State,I,Γ0,X,B0,BX,Par = Workspace
    @unpack  NUnique = Par


    Tol_Vertex = 1E16 #Initial tolerance level

    iG(x,nw) = iG_(State.γ,x,Lam,nw,Par)

    function getProp!(BubbleProp,nw1,nw2)
        for i in 1:NUnique, j in 1:NUnique
            BubbleProp[i,j] = iG(i,nw1) *iG(j,nw2)* Par.T
        end
        return BubbleProp
    end

    iter = 0
    while Tol_Vertex > Par.accuracy 
        iter+=1
        
        if iter > maxiterBubble
            @warn("BSE: No convergence found after $maxiterBubble iterations")
            break
        end

        writeTo!(OldState.Γ,State.Γ)

        setZero!( (B0,BX) )
                
        computeLeft2PartBubble!(B0,Γ0,Γ0,State.Γ,getProp!,Par,PropsBuffers,VertexBuffers)
        computeLeft2PartBubble!(BX,X,X,State.Γ,getProp!,Par,PropsBuffers,VertexBuffers)
        symmetrizeBubble!(B0,Par)
        symmetrizeBubble!(BX,Par)

        getXFromBubbles!(X,B0,BX) #TODO when applicable, this needs to be generalized for beyond-Parquet approximations 
        getVertexFromChannels!(State.Γ,I,X)
        
        Tol_Vertex = dist(OldState.Γ,State.Γ)
        
        iterateSDE!(Workspace,Lam,maxitergamma)
        CurrentObs = getObsFunc(Workspace,Lam)
        push!(Obs,CurrentObs)

        println("iteration $iter:")
        writeOutput(State,CurrentObs,Lam,Par)
        println("""
        Tol_Vertex = $Tol_Vertex
        """)
        isnan(Tol_Vertex) && @warn ": BSE: Solution diverged after $iter iterations\n"
    end
    return Workspace,Obs
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


function iterateSDE!(Workspace,Lam, maxiter)
    @unpack OldState,State,I,X,B0,BX,Par = Workspace
    Prop(x,nw) = -1/6*iG_(OldState.γ,x,Lam,nw,Par)

    Tol_gamma = 1E16
    iter = 0
    while Tol_gamma > 1E-5* Par.accuracy
        if iter >= maxiter
            # @warn("SDE: No convergence found after $maxiter iterations\nremaining Tol: $Tol_gamma")
            break
        end
        iter +=1
        writeTo!(OldState.γ,State.γ)
        setZero!(State.γ)
        compute1PartBubble!(State.γ,B0,Prop,Par)
        
        # @views println(State.γ[1:5])
        # @views println( Tuple(Prop(1,i) for i in 1:5 ))
        Tol_gamma = dist(State.γ,OldState.γ)
    end
    println("""
    \t\tSDE step done after $iter / $maxiter iterations """)
    return 
end

function getObservables(Workspace,Lam)
    @unpack MinimalOutput,N,np_vec,T,usesymmetry = Workspace.Par
    State = Workspace.State
    chi = getChi(Workspace.State,Lam,Workspace.Par)
    MaxVa = @view maximum(abs,State.Γ.a,dims = (2,3,4,5))[:,1,1,1]
    MaxVb = @view maximum(abs,State.Γ.b,dims = (2,3,4,5))[:,1,1,1]
    MaxVc = @view maximum(abs,State.Γ.c,dims = (2,3,4,5))[:,1,1,1]
    return Observables(chi,copy(State.γ),copy(State.f_int),MaxVa,MaxVb,MaxVc) # make sure to allocate new memory each time this function is called
end

function getChi(State::StateType, Lam::double,Par::Params)
	@unpack T,N,Npairs,lenIntw_acc,np_vec,invpairs,PairTypes,OnsitePairs = Par
	gamma = State.γ
	Vc = State.Γ.c

	iG(x,w) = iG_(gamma,x, Lam,w,Par)
	Vc_(Rij,s,t,u) = V_(Vc,Rij,s,t,u,invpairs[Rij],N)

	Chi = zeros(Npairs)

	@inbounds Threads.@threads for Rij in 1:Npairs
		@unpack xi,xj = PairTypes[Rij]
		for nK in -lenIntw_acc:lenIntw_acc-1
			if Rij in OnsitePairs
				Chi[Rij,1] += T * iG(xi,nK) ^2
			end
			for nK2 in -lenIntw_acc:lenIntw_acc-1
				npwpw2 = nK+nK2+1
				wmw2 = nK-nK2
				#use that Vc_0 is calculated from Vb
				GGGG = iG(xi,nK)^2 * iG(xj,nK2)^2
				Chi[Rij] += T^2 * GGGG *Vc_(Rij,0,npwpw2,wmw2)
			end
        end
    end
	return(Chi)
end


# writeOutput(St::StateType,Obs,Lam,Par) = println(Obs)
writeOutput(St::StateType,Obs,Lam,Par) = writeOutput(St.f_int,St.γ,St.Γ.a,St.Γ.b,St.Γ.c,Obs,Lam,Par)