
"""
Computes a two-particle bubble in the s-Channel given two four-point functions (i.e. vertices or other bubbles).
Γ is assumed to be a vertex.
To allow the computation using just the left part of a bubble, specification of the transpose is needed i.e. Transpose(X_L) = X_R and Transpose(X) = X, where X = XL + XR is the full bubble.
"""
function computeLeft2PartBubble!(ResultBubble,X,XTransp,Γ,getProp!,Par,Buffer)
    @unpack T,N,lenIntw,np_vec = Par.NumericalParams
    setZero!(ResultBubble)
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
                        addBLTilde!(ResultBubble,X,XTransp,Γ, is,it,iu,nw,Par,sprop)
                        if(!Par.Options.usesymmetry || nu<=nt)
                            addBL!(ResultBubble,X,XTransp,Γ,is,it,iu,nw,Par,sprop,VBuffer,XBuffer)
                        end
					end
				end
			end
		end
	end
    symmetrizeBubble!(ResultBubble,Par)
    return ResultBubble
end

"""
Computes a single-particle (i.e. self-energy) bubble. Allows specification of function type, i.e. what vertices are used since this is different if a bubble function is inserted as opposed to a vertex.
"""
function _compute1PartBubble!(Dgamma,XT1_::Function,XT2_::Function,Prop,Par)

    @unpack T,N,Ngamma,lenIntw_acc,np_vec_gamma = Par.NumericalParams
    @unpack siteSum,invpairs,Nsum,OnsitePairs = Par.System

    setZero!(Dgamma)
	Threads.@threads for iw1 in 1:Ngamma
		nw1 = np_vec_gamma[iw1]
    	for (x,Rx) in enumerate(OnsitePairs)
			for nw in -lenIntw_acc: lenIntw_acc-1
				jsum = 0.
				w1pw = nw1+nw+1 #w1 + w: Adding two Matsubara frequencies gives a +1
				w1mw = nw1-nw
				for k_spl in 1:Nsum[Rx]
					@unpack m,ki,xk = siteSum[k_spl,Rx]
					ik = invpairs[ki] # pair ik is inversed relative to pre-generated site summation in X (ki)!
					jsum += (XT1_(ik,w1pw,0,w1mw)+2*XT2_(ik,w1pw,0,w1mw))*Prop(xk,nw)*m
				end
				Dgamma[x,iw1] += T *jsum #For the self-energy derivative, the factor of 1/2 must be in the propagator
			end
		end
	end
    return Dgamma
end

"""
Computes a single-particle (i.e. self-energy) bubble. Can only be used if B is a bubble function
"""
function compute1PartBubble!(Dgamma,X::BubbleType,Prop,Par)
	XTa_(Rij,s,t,u) = XT_(X.a,X.a,Rij,s,t,u,Par.System.invpairs[Rij],Par.NumericalParams.N)
	XTc_(Rij,s,t,u) = XT_(X.c,X.b,Rij,s,t,u,Par.System.invpairs[Rij],Par.NumericalParams.N)
    _compute1PartBubble!(Dgamma,XTa_,XTc_,Prop,Par)
end

"""
Computes a single-particle (i.e. self-energy) bubble. Can only be used if argument is a vertex
"""
function compute1PartBubble!(Dgamma,Γ::VertexType,Prop,Par)
    @warn "compute1PartBubble! for vertices is not tested yet!"
	ΓTa_(Rij,s,t,u) = V_(Γ.a,Rij,t,u,s,Par.System.invpairs[Rij],Par.NumericalParams.N) # Tilde-type can be obtained by permutation of vertices
	ΓTc_(Rij,s,t,u) = V_(Γ.b,Rij,t,u,s,Par.System.invpairs[Rij],Par.NumericalParams.N) # cTilde corresponds to b type vertex!
    _compute1PartBubble!(Dgamma,Γa_,Γb_,Prop,Par)
end

@inline function X_(X::AbstractArray,XTransp::AbstractArray, Rj::Integer, ns::Integer,nt::Integer,nu::Integer,Rji::Integer,N::Integer)
    # @assert (ns+nt+nu) %2 != 0 "$ns + $nt +  $nu = $(ns+nt+nu)"
    ns,nt,nu,swapsites = convertFreqArgs(ns,nt,nu,N)
    Vertex = ifelse(swapsites,XTransp,X)
    Rj = ifelse(swapsites,Rji,Rj)
    return @inbounds Vertex[Rj,ns+1,nt+1,nu+1]
end

@inline function XT_(XT::AbstractArray,XTTransp::AbstractArray, Rj::Integer, ns::Integer,nt::Integer,nu::Integer,Rji::Integer,N::Integer)
    # @assert (ns+nt+nu) %2 != 0 "$ns + $nt +  $nu = $(ns+nt+nu)"
    Vertex = ifelse(nt*nu<0,XTTransp,XT)
    ns,nt,nu,swapsites = convertFreqArgsXT(ns,nt,nu,N)
    Rj = ifelse(swapsites,Rji,Rj)
    return @inbounds Vertex[Rj,ns+1,nt+1,nu+1]
end

@inline function bufferXT_!(Cache, X::AbstractArray, XTransp::AbstractArray, ns::Integer,nt::Integer,nu::Integer,invpairs::AbstractArray,N)
    
    Vertex = ifelse(nt*nu<0,XTransp,X)
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

@inline function fillBuffer!(VBuffer::VertexBufferType,XBuffer::BubbleBufferType,XL::BubbleType,XR::BubbleType,Γ::VertexType,is::Integer, it::Integer, iu::Integer, nwpr::Integer,Par::PMFRGParams)
    invpairs = Par.System.invpairs
    @unpack N,np_vec = Par.NumericalParams

    @unpack Va34,Vb34,Vc34,Vc43 = VBuffer
    @unpack XTa21,XTb21,XTc21,XTd21 = XBuffer

    ns = np_vec[is]
	nt = np_vec[it]
	nu = np_vec[iu]
	wpw1,wpw2,wmw3,wmw4 = mixedFrequencies(ns,nt,nu,nwpr)

	bufferV_!(Va34, Γ.a, ns, wmw3, wmw4, invpairs, N)
	bufferV_!(Vb34, Γ.b, ns, wmw3, wmw4, invpairs, N)
	bufferV_!(Vc34, Γ.c, ns, wmw3, wmw4, invpairs, N)
	
	bufferV_!(Vc43, Γ.c, ns, wmw4, wmw3, invpairs, N)

    bufferXT_!(XTa21, XR.Ta, XL.Ta, wpw2, ns, wpw1, invpairs, N)
    bufferXT_!(XTb21, XR.Tb, XL.Tb, wpw2, ns, wpw1, invpairs, N)
    bufferXT_!(XTc21, XR.Tc, XL.Tc, wpw2, ns, wpw1, invpairs, N)
    bufferXT_!(XTd21, XR.Td, XL.Td, wpw2, ns, wpw1, invpairs, N)
end

@inline function fillBuffer!(VBuffer::VertexBufferType,Γ0::BareVertexType,Γ::VertexType,is::Integer, it::Integer, iu::Integer, nwpr::Integer,Par::PMFRGParams)
    
    invpairs = Par.System.invpairs
    @unpack N,np_vec = Par.NumericalParams

    @unpack Va34,Vb34,Vc34,Vc43 = VBuffer
    ns = np_vec[is]
	nt = np_vec[it]
	nu = np_vec[iu]
	wpw1,wpw2,wmw3,wmw4 = mixedFrequencies(ns,nt,nu,nwpr)

	bufferV_!(Va34, Γ.a, ns, wmw3, wmw4, invpairs, N)
	bufferV_!(Vb34, Γ.b, ns, wmw3, wmw4, invpairs, N)
	bufferV_!(Vc34, Γ.c, ns, wmw3, wmw4, invpairs, N)
	
	bufferV_!(Vc43, Γ.c, ns, wmw4, wmw3, invpairs, N)
end

"""
adds to ResultBubble given the vertex as well as a bubble inserted on the left. Assumes that vertices and bubbles are given pre-computed in VertexBuffer.
"""
@inline function addBL!(B::BubbleType,XL::BubbleType,XR::BubbleType,Γ::VertexType,is::Integer, it::Integer, iu::Integer, nwpr::Integer,Par::PMFRGParams,Props,VBuffer::VertexBufferType,XBuffer::BubbleBufferType)
    @unpack N,np_vec = Par.NumericalParams
    @unpack Npairs,Nsum,siteSum,invpairs = Par.System

    fillBuffer!(VBuffer,XBuffer,XL,XR,Γ,is,it,iu,nwpr,Par)

    @unpack Va34,Vb34,Vc34,Vc43 = VBuffer
    @unpack XTa21,XTb21,XTc21,XTd21 = XBuffer

    S_ki = siteSum.ki
	S_kj = siteSum.kj
	S_xk = siteSum.xk
	S_m = siteSum.m

    @inbounds for Rij in 1:Npairs
        #loop over all left hand side inequivalent pairs Rij
        Ba_sum = 0. #Perform summation on this temp variable before writing to State array as Base.setindex! proved to be a bottleneck!
        Bb_sum = 0.
        Bc_sum = 0.
        @turbo unroll = 1 for k_spl in 1:Nsum[Rij]
            #loop over all Nsum summation elements defined in geometry. This inner loop is responsible for most of the computational effort! 
            ki,kj,m,xk = S_ki[k_spl,Rij],S_kj[k_spl,Rij],S_m[k_spl,Rij],S_xk[k_spl,Rij]

            Ptm = Props[xk,xk]*m
            
            Ba_sum += (
                Va34[kj] * XTa21[ki] + 
                2*Vb34[kj] * XTc21[ki] 
            )* Ptm

            Bb_sum += (
                Vb34[kj] * XTa21[ki] + 
                Va34[kj] * XTc21[ki] + 
                Vb34[kj] * XTc21[ki]
            )* Ptm
        end
        #split this into two loops because @turbo randomly allocates memory if expression is too long :(
        @turbo unroll = 1 for k_spl in 1:Nsum[Rij]
            #loop over all Nsum summation elements defined in geometry. This inner loop is responsible for most of the computational effort! 
            ki,kj,m,xk = S_ki[k_spl,Rij],S_kj[k_spl,Rij],S_m[k_spl,Rij],S_xk[k_spl,Rij]

            Ptm = Props[xk,xk]*m

            Bc_sum += (
                -Vc43[kj] * XTb21[ki] - 
                Vc34[kj] * XTd21[ki]
            )* Ptm
        end
        B.a[Rij,is,it,iu] += Ba_sum
        B.b[Rij,is,it,iu] += Bb_sum
        B.c[Rij,is,it,iu] += Bc_sum
    end
    return
end

@inline function addBL!(B::BubbleType,Γ0::BareVertexType,Γ::VertexType,is::Integer, it::Integer, iu::Integer, nwpr::Integer,Par::PMFRGParams,Props,Buffer::VertexBufferType)
    @unpack N,np_vec = Par.NumericalParams
    @unpack Npairs,Nsum,siteSum,invpairs = Par.System

    fillBuffer!(Buffer,Γ0,Γ,is,it,iu,nwpr,Par)
    
    @unpack Va34,Vb34,Vc34,Vc43 = Buffer
    S_ki = siteSum.ki
	S_kj = siteSum.kj
	S_xk = siteSum.xk
	S_m = siteSum.m

    @inbounds for Rij in 1:Npairs
        Bc_sum = 0.
        @turbo unroll = 1 for k_spl in 1:Nsum[Rij]
            #loop over all Nsum summation elements defined in geometry. This inner loop is responsible for most of the computational effort! 
            ki,kj,m,xk = S_ki[k_spl,Rij],S_kj[k_spl,Rij],S_m[k_spl,Rij],S_xk[k_spl,Rij]

            Ptm = Props[xk,xk]*m

            Bc_sum += (
                Vc43[kj] * Γ0.c[ki]+
                Vc34[kj] * Γ0.c[ki]
            )* Ptm
        end
        B.c[Rij,is,it,iu] += Bc_sum
    end
    return
end

@inline addBL!(B::BubbleType,Γ0L::BareVertexType,Γ0R::BareVertexType,Γ::VertexType,is::Integer, it::Integer, iu::Integer, nwpr::Integer,Par::PMFRGParams,Props,Buffer::VertexBufferType,XB::BubbleBufferType) = addBL!(B,Γ0L,Γ,is, it, iu, nwpr,Par,Props,Buffer)
##

function addBLTilde!(B::BubbleType,XL::BubbleType,XR::BubbleType,Γ::VertexType, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Par::PMFRGParams,Props)
    
    @unpack Npairs,invpairs,PairTypes,OnsitePairs = Par.System
    @unpack N,np_vec = Par.NumericalParams

    @inline Va_(Rij,s,t,u) = V_(Γ.a,Rij,s,t,u,invpairs[Rij],N)
    @inline Vb_(Rij,s,t,u) = V_(Γ.b,Rij,s,t,u,invpairs[Rij],N)
    @inline Vc_(Rij,s,t,u) = V_(Γ.c,Rij,s,t,u,invpairs[Rij],N)

    @inline XLa_(Rij,s,t,u) = X_(XL.a,XR.a,Rij,s,t,u,invpairs[Rij],N)
    @inline XLb_(Rij,s,t,u) = X_(XL.b,XR.b,Rij,s,t,u,invpairs[Rij],N)
    @inline XLc_(Rij,s,t,u) = X_(XL.c,XR.c,Rij,s,t,u,invpairs[Rij],N)

    @inline XRTa_(Rij,s,t,u) = XT_(XR.Ta,XL.Ta,Rij,s,t,u,invpairs[Rij],N)
    @inline XRTb_(Rij,s,t,u) = XT_(XR.Tb,XL.Tb,Rij,s,t,u,invpairs[Rij],N)
    @inline XRTc_(Rij,s,t,u) = XT_(XR.Tc,XL.Tc,Rij,s,t,u,invpairs[Rij],N)
    @inline XRTd_(Rij,s,t,u) = XT_(XR.Td,XL.Td,Rij,s,t,u,invpairs[Rij],N)

	ns = np_vec[is]
	nt = np_vec[it]
	nu = np_vec[iu]
	wpw1,wpw2,wmw3,wmw4 = mixedFrequencies(ns,nt,nu,nwpr)

    #Btilde only defined for nonlocal pairs Rij != Rii
    for Rij in 1:Npairs
        Rij in OnsitePairs && continue
        #loop over all left hand side inequivalent pairs Rij
        Rji = invpairs[Rij] # store pair corresponding to Rji (easiest case: Rji = Rij)
        @unpack xi,xj = PairTypes[Rij]

        B.Ta[Rij,is,it,iu] += 
        Props[xj, xi]*(
        Va_(Rij, wmw4, ns, wmw3)*XLa_(Rij, wpw2, ns, wpw1) + 
        2*Vc_(Rij, wmw4, ns, wmw3)*XLc_(Rij, wpw2, ns, wpw1)) + 
        Props[xi, xj]*(
        Va_(Rji, wmw3, ns, wmw4)*XRTa_(Rji, wpw2, wpw1, ns) + 
        2*Vc_(Rji, wmw3, ns, wmw4)*XRTd_(Rji, wpw2, wpw1, ns))
        
        B.Tb[Rij,is,it,iu] += 
        Props[xj, xi]*(
        Va_(Rij, wmw4, ns, wmw3)*XLc_(Rij, wpw2, ns, wpw1) + 
        Vc_(Rij, wmw4, ns, wmw3)*(
        XLa_(Rij, wpw2, ns, wpw1) + XLc_(Rij, wpw2, ns, wpw1))) + 
        Props[xi, xj]*(
        Va_(Rji, wmw3, ns, wmw4)*XRTd_(Rji, wpw2, wpw1, ns) + 
        Vc_(Rji, wmw3, ns, wmw4)*(
        XRTa_(Rji, wpw2, wpw1, ns) + XRTd_(Rji, wpw2, wpw1, ns)))

        B.Tc[Rij,is,it,iu] += 
        Props[xj, xi]*(
        Vb_(Rij, wmw4, ns, wmw3)*XLb_(Rij, wpw2, ns, wpw1) + 
        Vc_(Rij, wmw4, wmw3, ns)*XLc_(Rij, wpw2, wpw1, ns)) + 
        Props[xi, xj]*(
        -Vc_(Rji, wmw3, wmw4, ns)*XRTb_(Rji, wpw2, wpw1, ns) + 
        Vb_(Rji, wmw3, ns, wmw4)*XRTc_(Rji, wpw2, wpw1, ns))
    end
end
##

function addBLTilde!(B::BubbleType,Γ0::BareVertexType,Γ::VertexType, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Par::PMFRGParams,Props)

    @unpack Npairs,invpairs,PairTypes,OnsitePairs = Par.System
    @unpack N,np_vec = Par.NumericalParams


    @inline Va_(Rij,s,t,u) = V_(Γ.a,Rij,s,t,u,invpairs[Rij],N)
    @inline Vb_(Rij,s,t,u) = V_(Γ.b,Rij,s,t,u,invpairs[Rij],N)
    @inline Vc_(Rij,s,t,u) = V_(Γ.c,Rij,s,t,u,invpairs[Rij],N)

	ns = np_vec[is]
	nt = np_vec[it]
	nu = np_vec[iu]
	wpw1,wpw2,wmw3,wmw4 = mixedFrequencies(ns,nt,nu,nwpr)

    #Btilde only defined for nonlocal pairs Rij != Rii
    for Rij in 1:Npairs
        Rij in OnsitePairs && continue
        #loop over all left hand side inequivalent pairs Rij
        Rji = invpairs[Rij] # store pair corresponding to Rji (easiest case: Rji = Rij)
        @unpack xi,xj = PairTypes[Rij]

        B.Ta[Rij,is,it,iu] += 
        Props[xj, xi]*2*Vc_(Rij, wmw4, ns, wmw3)*Γ0.c[Rij] + 
        Props[xi, xj]*2*Vc_(Rji, wmw3, ns, wmw4)*Γ0.c[Rji]
        
        B.Tb[Rij,is,it,iu] += 
        Props[xj, xi]*(Va_(Rij, wmw4, ns, wmw3) + Vc_(Rij, wmw4, ns, wmw3))*Γ0.c[Rij]
        Props[xi, xj]*(Va_(Rji, wmw3, ns, wmw4) + Vc_(Rji, wmw3, ns, wmw4))*Γ0.c[Rji]

        B.Tc[Rij,is,it,iu] += 
        Props[xj, xi]*Vc_(Rij, wmw4, wmw3, ns)*Γ0.c[Rij] +
        Props[xi, xj]*Vc_(Rji, wmw3, ns, wmw4)*Γ0.c[Rji]
    end
end
@inline addBLTilde!(B::BubbleType,Γ0L::BareVertexType,Γ0R::BareVertexType,Γ::VertexType, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Par::PMFRGParams,Props) = addBLTilde!(B,Γ0L,Γ, is, it, iu, nwpr, Par, Props)