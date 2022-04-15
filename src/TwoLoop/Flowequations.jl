function getDeriv!(Deriv,State,setup::Tuple{BubbleType,BubbleType,T,TwoLoopParams},Lam) where T
        
    X,Y,Buffs,Par = setup #use pre-allocated X and XTilde to reduce garbage collector time
    Workspace = TwoLoopWorkspace(Deriv,State,X,Y,Buffs,Par)

    getDFint!(Workspace,Lam)
    get_Self_Energy!(Workspace,Lam)
    getXBubble!(Workspace,Lam)
    symmetrizeBubble!(Workspace.X,Par)
    
    addToVertexFromBubble!(Workspace.Deriv.Γ,Workspace.X)
    getTwoLoopDeriv!(Workspace,Lam)
    
    symmetrizeBubble!(Workspace.Y,Par)
    addToVertexFromBubble!(Workspace.Deriv.Γ,Workspace.Y)

    symmetrizeVertex!(Workspace.Deriv.Γ,Par)
    return
end

function getDerivVerbose!(Deriv,State,setup::Tuple{BubbleType,BubbleType,T,TwoLoopParams},Lam)
    X,Y,Buffs,Par = setup #use pre-allocated X and XTilde to reduce garbage collector time
    print("Workspace:\n\t") 
    @time Workspace = TwoLoopWorkspace(Deriv,State,X,Y,Buffs,Par)

    print("getDFint:\n\t") 
    @time getDFint!(Workspace,Lam)

    print("get_Self_Energy:\n\t") 
    @time get_Self_Energy!(Workspace,Lam)

    print("getVertexDeriv:\n\t") 
    @time getXBubble!(Workspace,Lam)

    print("SymmetryX:\n\t") 
    symmetrizeBubble!(Workspace.X,Par)

    print("TwoLoop:\n\t") 
    @time getTwoLoopDeriv!(Workspace,Lam)
    
    print("SymmetryY:\n\t") 
    @time symmetrizeBubble!(Workspace.Y,Par)
    
    symmetrizeVertex!(Workspace.State.Γ,Par)
    return
end

function getTwoLoopDeriv!(Workspace::TwoLoopWorkspace,Lam)
    @unpack X,Y,State,Par,Buffer = Workspace
    getProp! = constructPropagatorFunction(Workspace,Lam)
    @unpack State,Par = Workspace
    @unpack T,N,np_vec,lenIntw = Par.NumericalParams
    @unpack NUnique,OnsitePairs = Par.System
    setZero!(Y)
    addToLeft2PartBubble!(Y,X,X,State.Γ,getProp!,Par,Buffer)
    addToRight2PartBubble!(Y,X,X,State.Γ,getProp!,Par,Buffer)
end

@inline function convertFreqArgsXT(ns,nt,nu,Nw)
    # @assert (ns+nt+nu) %2 != 0 "trying to convert wrong freqs $ns + $nt +  $nu = $(ns+nt+nu)"
    swapsites = ns*nu <0
    ns,nt,nu = abs.((ns,nt,nu))
    ns = min( ns, Nw - 1 - (ns+Nw-1)%2)
    nt = min( nt, Nw - 1 - (nt+Nw-1)%2)
    nu = min( nu, Nw - 1 - (nu+Nw-1)%2)
    return ns,nt,nu,swapsites
end

"""
Computes a two-particle bubble in the s-Channel given two four-point functions (i.e. vertices or other bubbles).
Γ is assumed to be a vertex.
To allow the computation using just the left part of a bubble, specification of the transpose is needed i.e. Transpose(X_L) = X_R and Transpose(X) = X, where X = XL + XR is the full bubble.
"""
function addTo2PartBubble!(ResultBubble::BubbleType,X,XTransp,Γ::VertexType,getProp!::Function,addBTilde_Func!::Function,addB_Func!::Function,Par,Buffer)
    @unpack T,N,lenIntw,np_vec = Par.NumericalParams
    @sync begin
		for is in 1:N,it in 1:N
			Threads.@spawn begin
				BubbleProp = Buffer.Props[Threads.threadid()] # get pre-allocated thread-safe buffers
				VBuffer = Buffer.Vertex[Threads.threadid()]
				XBuffer = Buffer.X[Threads.threadid()]
				ns = np_vec[is]
				nt = np_vec[it]
				for iu in 1:N
					nu = np_vec[iu]
					if (ns+nt+nu)%2 == 0	# skip unphysical bosonic frequency combinations
						continue
					end
					for nw in -lenIntw:lenIntw-1 # Matsubara sum
                        sprop = getProp!(BubbleProp,nw,nw+ns)
                        addBTilde_Func!(ResultBubble,X,XTransp,Γ, is,it,iu,nw,Par,sprop)
                        if(!Par.Options.usesymmetry || nu<=nt)
                            addB_Func!(ResultBubble,X,XTransp,Γ,is,it,iu,nw,Par,sprop,VBuffer,XBuffer)
                        end
					end
				end
			end
		end
	end
    return ResultBubble
end

@inline function compute2PartBubble!(ResultBubble::BubbleType,X,XTransp,Γ::VertexType,getProp!::Function,addBTilde_Func!::Function,addB_Func!::Function,Par,Buffer)
    setZero!(ResultBubble)
    addTo2PartBubble!(ResultBubble,X,XTransp,Γ,getProp!,addBTilde_Func!,addB_Func!,Par,Buffer)
    symmetrizeBubble!(ResultBubble,Par)
end

@inline addToLeft2PartBubble!(ResultBubble::BubbleType,X,XTransp,Γ::VertexType,getProp!::Function,Par,Buffer) = addTo2PartBubble!(ResultBubble,X,XTransp,Γ,getProp!,addBLTilde!,addBL!,Par,Buffer)

@inline addToRight2PartBubble!(ResultBubble::BubbleType,X,XTransp,Γ::VertexType,getProp!::Function,Par,Buffer) = addTo2PartBubble!(ResultBubble,X,XTransp,Γ,getProp!,addBRTilde!,addBR!,Par,Buffer)

@inline computeLeft2PartBubble!(ResultBubble::BubbleType,X,XTransp,Γ::VertexType,getProp!::Function,Par,Buffer) = compute2PartBubble!(ResultBubble,X,XTransp,Γ,getProp!,addBLTilde!,addBL!,Par,Buffer)

@inline computeRight2PartBubble!(ResultBubble::BubbleType,X,XTransp,Γ::VertexType,getProp!::Function,Par,Buffer) = compute2PartBubble!(ResultBubble,X,XTransp,Γ,getProp!,addBRTilde!,addBR!,Par,Buffer)