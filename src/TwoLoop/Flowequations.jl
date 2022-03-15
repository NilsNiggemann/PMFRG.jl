function getDeriv2L!(Deriv,State,setup,Lam)
        
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

function getDeriv2LVerbose!(Deriv,State,setup,Lam)
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
    
    print("SymmetryX:\n\t") 
    @time symmetrizeBubble!(Workspace.Y,Par)
    
    symmetrizeVertex!(Workspace.State.Γ,Par)
    return
end

function getTwoLoopDeriv!(Workspace::TwoLoopWorkspace,Lam)
    @unpack State,Par = Workspace
    @unpack T,N,np_vec,lenIntw = Par.NumericalParams
    @unpack NUnique,OnsitePairs = Par.System
    
    PropsBuffers = Workspace.Buffer.Props # Todo make Buffer a Structarray?
    VertexBuffers = Workspace.Buffer.Vertex
    XBuffers = Workspace.Buffer.X

    iG(x,nw) = iG_(State.γ,x,Lam,nw,T)

    function getProp!(BubbleProp,nw1,nw2)
        for i in 1:NUnique, j in 1:NUnique
            BubbleProp[i,j] = iG(i,nw1) *iG(j,nw2)* T
        end
        return BubbleProp
    end
    @sync begin
		for is in 1:N,it in 1:N
			Threads.@spawn begin
				BubbleProp = PropsBuffers[Threads.threadid()] # get pre-allocated thread-safe buffers
				VBuffer = VertexBuffers[Threads.threadid()]
				XBuffer = XBuffers[Threads.threadid()]
				ns = np_vec[is]
				nt = np_vec[it]
				for iu in 1:N
					nu = np_vec[iu]
					if (ns+nt+nu)%2 == 0	# skip unphysical bosonic frequency combinations
						continue
					end
					for nw in -lenIntw:lenIntw-1 # Matsubara sum
                        sprop = getProp!(BubbleProp,nw,nw+ns)
                        addYTilde!(Workspace,is,it,iu,nw,sprop)
                        if(!Par.Options.usesymmetry || nu<=nt)
                            addY!(Workspace,is,it,iu,nw,sprop,VBuffer,XBuffer)
                        end
					end
				end
			end
		end
	end
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

function XT_(Vertex::AbstractArray, Rj::Integer, ns::Integer,nt::Integer,nu::Integer,Rji::Integer,N::Integer)
    # @assert (ns+nt+nu) %2 != 0 "$ns + $nt +  $nu = $(ns+nt+nu)"
    ns,nt,nu,swapsites = convertFreqArgsXT(ns,nt,nu,N)
    Rj = ifelse(swapsites,Rji,Rj)
    return @inbounds Vertex[Rj,ns+1,nt+1,nu+1]
end

@inline function bufferXT_!(Cache, Vertex::AbstractArray, ns::Integer,nt::Integer,nu::Integer,invpairs::AbstractArray,N)
    ns,nt,nu,swapsites = convertFreqArgsXT(ns,nt,nu,N)
    # @assert (ns+nt+nu) %2 != 0 "$ns + $nt +  $nu = $(ns+nt+nu)"

    is,it,iu = ns+1,nt+1,nu+1
    @inbounds begin 
        if swapsites
            @turbo unroll = 1 inline = true for R in eachindex(Cache,invpairs)
                Cache[R] = Vertex[invpairs[R],is,it,iu]
            end
        else
            @turbo unroll = 1 inline = true for R in eachindex(Cache,invpairs)
                Cache[R] = Vertex[R,is,it,iu]
            end
        end
    end
end

"""
like addX! but with Y
"""
@inline function addY!(Workspace::PMFRGWorkspace, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Props,VertexBuffer,XBuffer)
	@unpack State,X,Y,Par = Workspace
	@unpack N,np_vec = Par.NumericalParams
	@unpack Npairs,Nsum,siteSum,invpairs = Par.System
    
    @unpack Va12,Vb12,Vc12,Va34,Vb34,Vc34,Vc21,Vc43 = VertexBuffer
    @unpack XTa21,XTa43,XTb21,XTb43,XTc21,XTc43,XTd21,XTd43 = XBuffer

    ns = np_vec[is]
	nt = np_vec[it]
	nu = np_vec[iu]
	wpw1,wpw2,wmw3,wmw4 = mixedFrequencies(ns,nt,nu,nwpr)

    bufferV_!(Va12, State.Γ.a , ns, wpw1, wpw2, invpairs, N)
	bufferV_!(Vb12, State.Γ.b , ns, wpw1, wpw2, invpairs, N)
	bufferV_!(Vc12, State.Γ.c , ns, wpw1, wpw2, invpairs, N)

	bufferV_!(Va34, State.Γ.a , ns, wmw3, wmw4, invpairs, N)
	bufferV_!(Vb34, State.Γ.b , ns, wmw3, wmw4, invpairs, N)
	bufferV_!(Vc34, State.Γ.c , ns, wmw3, wmw4, invpairs, N)
	
	bufferV_!(Vc21, State.Γ.c , ns, wpw2, wpw1, invpairs, N)
	bufferV_!(Vc43, State.Γ.c , ns, wmw4, wmw3, invpairs, N)

    bufferXT_!(XTa21,X.Ta , wpw2, ns, wpw1, invpairs, N)
	bufferXT_!(XTa43,X.Ta , wmw4, ns, wmw3, invpairs, N)

    bufferXT_!(XTb21,X.Tb , wpw2, ns, wpw1, invpairs, N)
	bufferXT_!(XTb43,X.Tb , wmw4, ns, wmw3, invpairs, N)

    bufferXT_!(XTc21,X.Tc , wpw2, ns, wpw1, invpairs, N)
	bufferXT_!(XTc43,X.Tc , wmw4, ns, wmw3, invpairs, N)

    bufferXT_!(XTd21,X.Td , wpw2, ns, wpw1, invpairs, N)
	bufferXT_!(XTd43,X.Td , wmw4, ns, wmw3, invpairs, N)

    S_ki = siteSum.ki
	S_kj = siteSum.kj
	S_xk = siteSum.xk
	S_m = siteSum.m

    @inbounds for Rij in 1:Npairs
        #loop over all left hand side inequivalent pairs Rij
        Ya_sum = 0. #Perform summation on this temp variable before writing to State array as Base.setindex! proved to be a bottleneck!
        Yb_sum = 0.
        Yc_sum = 0.
        @turbo unroll = 1 for k_spl in 1:Nsum[Rij]
            #loop over all Nsum summation elements defined in geometry. This inner loop is responsible for most of the computational effort! 
            ki,kj,m,xk = S_ki[k_spl,Rij],S_kj[k_spl,Rij],S_m[k_spl,Rij],S_xk[k_spl,Rij]

            Ptm = Props[xk,xk]*m
            
            Ya_sum += (
                Va34[kj] * XTa21[ki] + 
                Va12[ki] * XTa43[kj] + 
                2*Vb34[kj] * XTc21[ki] + 
                2*Vb12[ki] * XTc43[kj]
            )* Ptm

            Yb_sum += (
                Vb34[kj] * XTa21[ki] + 
                Vb12[ki] * XTa43[kj] + 
                Va34[kj] * XTc21[ki] + 
                Vb34[kj] * XTc21[ki] + 
                Va12[ki] * XTc43[kj] + 
                Vb12[ki] * XTc43[kj]
            )* Ptm
        end
        #split this into two loops because @turbo randomly allocates memory if expression is too long :(
        @turbo unroll = 1 for k_spl in 1:Nsum[Rij]
            #loop over all Nsum summation elements defined in geometry. This inner loop is responsible for most of the computational effort! 
            ki,kj,m,xk = S_ki[k_spl,Rij],S_kj[k_spl,Rij],S_m[k_spl,Rij],S_xk[k_spl,Rij]

            Ptm = Props[xk,xk]*m

            Yc_sum += (
                -Vc43[kj] * XTb21[ki] - 
                Vc21[ki] * XTb43[kj] + 
                Vc34[kj] * XTd21[ki] + 
                Vc12[ki] * XTd43[kj]
            )* Ptm
        end
        Y.a[Rij,is,it,iu] += Ya_sum
        Y.b[Rij,is,it,iu] += Yb_sum
        Y.c[Rij,is,it,iu] += Yc_sum
    end
    return
end
##

function addYTilde!(Workspace::PMFRGWorkspace, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Props)
	@unpack State,X,Y,Par = Workspace
	@unpack N,np_vec = Par.NumericalParams
	@unpack Npairs,invpairs,PairTypes,OnsitePairs = Par.System


    @inline Va_(Rij,s,t,u) = V_(State.Γ.a,Rij,s,t,u,invpairs[Rij],N)
    @inline Vb_(Rij,s,t,u) = V_(State.Γ.b,Rij,s,t,u,invpairs[Rij],N)
    @inline Vc_(Rij,s,t,u) = V_(State.Γ.c,Rij,s,t,u,invpairs[Rij],N)

    @inline Xa_(Rij,s,t,u) = V_(X.a,Rij,s,t,u,invpairs[Rij],N)
    @inline Xb_(Rij,s,t,u) = V_(X.b,Rij,s,t,u,invpairs[Rij],N)
    @inline Xc_(Rij,s,t,u) = V_(X.c,Rij,s,t,u,invpairs[Rij],N)

    @inline XTa_(Rij,s,t,u) = XT_(X.Ta,Rij,s,t,u,invpairs[Rij],N)
    @inline XTb_(Rij,s,t,u) = XT_(X.Tb,Rij,s,t,u,invpairs[Rij],N)
    @inline XTc_(Rij,s,t,u) = XT_(X.Tc,Rij,s,t,u,invpairs[Rij],N)
    @inline XTd_(Rij,s,t,u) = XT_(X.Td,Rij,s,t,u,invpairs[Rij],N)

	ns = np_vec[is]
	nt = np_vec[it]
	nu = np_vec[iu]
	wpw1,wpw2,wmw3,wmw4 = mixedFrequencies(ns,nt,nu,nwpr)

    #Ytilde only defined for nonlocal pairs Rij != Rii
    for Rij in 1:Npairs
        Rij in OnsitePairs && continue
        #loop over all left hand side inequivalent pairs Rij
        Rji = invpairs[Rij] # store pair corresponding to Rji (easiest case: Rji = Rij)
        @unpack xi,xj = PairTypes[Rij]

        Y.Ta[Rij,is,it,iu] += 
        Props[xj, xi]*(
        Va_(Rij, wpw2, ns, wpw1)*Xa_(Rij, wmw4, ns, wmw3) + 
        Va_(Rij, wmw4, ns, wmw3)*Xa_(Rij, wpw2, ns, wpw1) + 
        2*Vc_(Rij, wpw2, ns, wpw1)*Xc_(Rij, wmw4, ns, wmw3) + 
        2*Vc_(Rij, wmw4, ns, wmw3)*Xc_(Rij, wpw2, ns, wpw1)) + 
        Props[xi, xj]*(
        Va_(Rji, wpw1, ns, wpw2)*XTa_(Rji, wmw4, wmw3, ns) + 
        Va_(Rji, wmw3, ns, wmw4)*XTa_(Rji, wpw2, wpw1, ns) + 
        2*Vc_(Rji, wpw1, ns, wpw2)*XTd_(Rji, wmw4, wmw3, ns) + 
        2*Vc_(Rji, wmw3, ns, wmw4)*XTd_(Rji, wpw2, wpw1, ns))
        
        Y.Tb[Rij,is,it,iu] += 
        Props[xj, xi]*(
        Va_(Rij, wpw2, ns, wpw1)*Xc_(Rij, wmw4, ns, wmw3) + 
        Vc_(Rij, wpw2, ns, wpw1)*(
        Xa_(Rij, wmw4, ns, wmw3) + Xc_(Rij, wmw4, ns, wmw3)) + 
        Va_(Rij, wmw4, ns, wmw3)*Xc_(Rij, wpw2, ns, wpw1) + 
        Vc_(Rij, wmw4, ns, wmw3)*(
        Xa_(Rij, wpw2, ns, wpw1) + Xc_(Rij, wpw2, ns, wpw1))) + 
        Props[xi, xj]*(
        Va_(Rji, wpw1, ns, wpw2)*XTd_(Rji, wmw4, wmw3, ns) + 
        Vc_(Rji, wpw1, ns, wpw2)*(
        XTa_(Rji, wmw4, wmw3, ns) + XTd_(Rji, wmw4, wmw3, ns)) + 
        Va_(Rji, wmw3, ns, wmw4)*XTd_(Rji, wpw2, wpw1, ns) + 
        Vc_(Rji, wmw3, ns, wmw4)*(
        XTa_(Rji, wpw2, wpw1, ns) + XTd_(Rji, wpw2, wpw1, ns)))

        Y.Tc[Rij,is,it,iu] += 
        Props[xj, xi]*(
        Vb_(Rij, wpw2, ns, wpw1)*Xb_(Rij, wmw4, ns, wmw3) + 
        Vb_(Rij, wmw4, ns, wmw3)*Xb_(Rij, wpw2, ns, wpw1) + 
        Vc_(Rij, wpw2, wpw1, ns)*Xc_(Rij, wmw4, wmw3, ns) + 
        Vc_(Rij, wmw4, wmw3, ns)*Xc_(Rij, wpw2, wpw1, ns)) + 
        Props[xi, xj]*(
        -Vc_(Rji, wpw1, wpw2, ns)*XTb_(Rji, wmw4, wmw3, ns) - 
        Vc_(Rji, wmw3, wmw4, ns)*XTb_(Rji, wpw2, wpw1, ns) + 
        Vb_(Rji, wpw1, ns, wpw2)*XTc_(Rji, wmw4, wmw3, ns) + 
        Vb_(Rji, wmw3, ns, wmw4)*XTc_(Rji, wpw2, wpw1, ns))
    end
end