struct ParquetWorkspace{T,P <: Params{S,M} where {S,M}}
    State::StateType{T} #Stores the current step of iteration
    NewState::StateType{T} #Stores the next step of iteration

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
    setToBareVertex!(Workspace.State.Γ,Par)
    return Workspace
end

"""Given set of parameters solve the self-consistent Parquet approximation iteratively. 
TODO: Save output to file"""
function SolveParquet(Par::Params,Lam;kwargs...)
    Workspace = SetupParquet(Par)
    
    PropsBuffers = [Matrix{double}(undef,Par.NUnique,Par.NUnique) for _ in 1:Threads.nthreads()] 
    VertexBuffers = [VertexBuffer(Par.Npairs) for _ in 1:Threads.nthreads()] 

    iterateSolution!(Workspace,Lam,PropsBuffers,VertexBuffers;kwargs...)
end

"""Obtains a solution to Bethe-Salpeter and Schwinger-Dyson equations by iteration until convergence is reached up to accuracy specified by accuracy in Params"""
function iterateSolution!(Workspace::ParquetWorkspace,Lam,PropsBuffers,VertexBuffers;maxiterBubble = 10, maxitergamma=100)
    @unpack State,NewState,I,Γ0,X,B0,BX,Par = Workspace
    @unpack  NUnique = Par


    Error_Vertex = 1E16 #Initial tolerance level

    iG(x,nw) = iG_(State.γ,x,Lam,nw,Par)

    function getProp!(BubbleProp,nw1,nw2)
        for i in 1:NUnique, j in 1:NUnique
            BubbleProp[i,j] = iG(i,nw1) *iG(j,nw2)* Par.T
        end
        return BubbleProp
    end

    iter = 0
    while Error_Vertex > Par.accuracy 
        iter+=1
        println("iteration $iter:")
        println("""
            Error_Vertex = $Error_Vertex
        """)
        if iter >= maxiterBubble
            @warn("No convergence of found after $iter iterations")
            break
        end

        setZero!( (B0,BX) )

        computeLeft2PartBubble!(B0,Γ0,Γ0,State.Γ,getProp!,Par,PropsBuffers,VertexBuffers)
        computeLeft2PartBubble!(BX,X,X,State.Γ,getProp!,Par,PropsBuffers,VertexBuffers)

        getXFromBubbles!(X,B0,BX) #TODO when applicable, this needs to be generalized for beyond-Parquet approximations 
        getVertexFromChannels!(NewState.Γ,I,X)

        Error_Vertex = dist(State.Γ,NewState.Γ)

        setEqual!(State.Γ,NewState.Γ)
        iterateSDE!(Workspace,Lam,maxiter = maxitergamma)
    end
    return Workspace
end

function getXFromBubbles!(X::BubbleType,B0::BubbleType,BX::BubbleType)
    for f in fieldnames(BubbleType)
        Xf,B0f,BXf = getfield(X,f),getfield(B0,f),getfield(BX,f)

        Threads.@threads for i in eachindex(Xf,B0f,BXf)
           Xf[i] = -0.5* B0f[i] + BXf[i]
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

function iterateSDE!(Workspace,Lam; maxiter = 100)
    @unpack State,NewState,I,X,B0,BX,Par = Workspace
    Prop(x,nw) = iG_(State.γ,x,Lam,nw,Par)

    Error_gamma = 1E16
    iter = 0
    while Error_gamma > Par.accuracy && iter < maxiter
        iter +=1
        if iter >= maxiter
            @warn("No convergence of SDE found after $iter iterations\n
            remaining error: $Error_gamma")
            break
        end
        compute1PartBubble!(NewState.γ,B0,Prop,Par)
        Error_gamma = dist(NewState.γ,State.γ)
        setEqual!(State.γ,NewState.γ)
    end
    println("""
    \t\t SDE step done after $iter iterations """)
    return 
end
