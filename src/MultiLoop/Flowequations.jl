

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

@inline function fillBufferL!(VBuffer::VertexBufferType,XBuffer::BubbleBufferType,XL::BubbleType,XR::BubbleType,Γ::VertexType,is::Integer, it::Integer, iu::Integer, nwpr::Integer,Par::PMFRGParams)
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


@inline function fillBufferR_new!(VBuffer::VertexBufferType,XBuffer::BubbleBufferType,XL::BubbleType,XR::BubbleType,Γ::VertexType,is::Integer, it::Integer, iu::Integer, nwpr::Integer,Par::PMFRGParams)
    invpairs = Par.System.invpairs
    @unpack N,np_vec = Par.NumericalParams

    @unpack Va34,Vb34,Vc34,Vc43 = VBuffer
    @unpack XTa21,XTb21,XTc21,XTd21 = XBuffer

    ns = np_vec[is]
	nt = np_vec[it]
	nu = np_vec[iu]
	wpw1,wpw2,wmw3,wmw4 = mixedFrequencies(ns,nt,nu,nwpr)
    wpw1,wpw2,wmw3,wmw4 = wmw3,wmw4,wpw1,wpw2 #swap frequencies for right bubble

	bufferV_!(Va34, Γ.a, ns, wmw3, wmw4, invpairs, N)
	bufferV_!(Vb34, Γ.b, ns, wmw3, wmw4, invpairs, N)
	bufferV_!(Vc34, Γ.c, ns, wmw3, wmw4, invpairs, N)
	
	bufferV_!(Vc43, Γ.c, ns, wmw4, wmw3, invpairs, N)

    bufferXT_!(XTa21, XR.Ta, XL.Ta, wpw2, ns, wpw1, invpairs, N)
    bufferXT_!(XTb21, XR.Tb, XL.Tb, wpw2, ns, wpw1, invpairs, N)
    bufferXT_!(XTc21, XR.Tc, XL.Tc, wpw2, ns, wpw1, invpairs, N)
    bufferXT_!(XTd21, XR.Td, XL.Td, wpw2, ns, wpw1, invpairs, N)
end

@inline function fillBufferR!(VBuffer::VertexBufferType,XBuffer::BubbleBufferType,XL::BubbleType,XR::BubbleType,Γ::VertexType,is::Integer, it::Integer, iu::Integer, nwpr::Integer,Par::PMFRGParams)
    invpairs = Par.System.invpairs
    @unpack N,np_vec = Par.NumericalParams

    @unpack Va12,Vb12,Vc12,Vc21 = VBuffer
    @unpack XTa43,XTb43,XTc43,XTd43 = XBuffer

    ns = np_vec[is]
	nt = np_vec[it]
	nu = np_vec[iu]
	wpw1,wpw2,wmw3,wmw4 = mixedFrequencies(ns,nt,nu,nwpr)


    bufferV_!(Va12, Γ.a , ns, wpw1, wpw2, invpairs, N)
	bufferV_!(Vb12, Γ.b , ns, wpw1, wpw2, invpairs, N)
	bufferV_!(Vc12, Γ.c , ns, wpw1, wpw2, invpairs, N)

	bufferV_!(Vc21, Γ.c , ns, wpw2, wpw1, invpairs, N)

	bufferXT_!(XTa43,XL.Ta,XR.Ta , wmw4, ns, wmw3, invpairs, N) # Caution check whether XR and XL do not need to be swapped
	bufferXT_!(XTb43,XL.Tb,XR.Tb , wmw4, ns, wmw3, invpairs, N) # Caution check whether XR and XL do not need to be swapped
	bufferXT_!(XTc43,XL.Tc,XR.Tc , wmw4, ns, wmw3, invpairs, N) # Caution check whether XR and XL do not need to be swapped
	bufferXT_!(XTd43,XL.Td,XR.Td , wmw4, ns, wmw3, invpairs, N) # Caution check whether XR and XL do not need to be swapped
end

@inline function fillBufferL!(VBuffer::VertexBufferType,Γ0::BareVertexType,Γ::VertexType,is::Integer, it::Integer, iu::Integer, nwpr::Integer,Par::PMFRGParams)
    
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
