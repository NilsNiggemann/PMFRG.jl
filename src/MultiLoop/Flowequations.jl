function getDeriv!(Deriv,State,setup,Lam)
        
    X,XTilde,Y,YTilde,PropsBuffers,VertexBuffers,Par = setup #use pre-allocated X and XTilde to reduce garbage collector time
    N = Par.N
    OneLoopWorkspace = Workspace_Struct(Deriv,State,X,XTilde)
    TwoLoopWorkspace = Y_Workspace_Struct(Y,YTilde)
    getDFint!(OneLoopWorkspace,Lam,Par)
    get_Self_Energy!(OneLoopWorkspace,Lam,Par)
    getVertexDeriv!(OneLoopWorkspace,Lam,Par,PropsBuffers,VertexBuffers)
    symmetrizeX!(OneLoopWorkspace,Par)
    getTwoLoopDeriv!(OneLoopWorkspace,TwoLoopWorkspace,Lam,Par,PropsBuffers,VertexBuffers)
    symmetrizeY!(OneLoopWorkspace,TwoLoopWorkspace,Par)
    return
end

function getDerivVerbose!(Deriv,State,XandPar,Lam)
    X,XTilde,Y,YTilde,Par = setup #use pre-allocated X and XTilde to reduce garbage collector time
    N = Par.N
    print("Workspace:\n\t") 
    @time begin
        OneLoopWorkspace = Workspace_Struct(Deriv,State,X,XTilde)
        TwoLoopWorkspace = Y_Workspace_Struct(Y,YTilde)
    end
    print("getDFint:\n\t") 
    @time getDFint!(Workspace,Lam,Par)

    print("get_Self_Energy:\n\t") 
    @time get_Self_Energy!(Workspace,Lam,Par)

    print("getVertexDeriv:\n\t") 
    @time getVertexDeriv!(Workspace,Lam,Par,PropsBuffers,VertexBuffers)

    print("SymmetryX:\n\t") 
    @time symmetrizeX!(Workspace,Par)

    print("TwoLoop:\n\t") 
    @time getTwoLoopDeriv!(OneLoopWorkspace,TwoLoopWorkspace,Lam,Par,PropsBuffers,VertexBuffers)
    
    print("SymmetryX:\n\t") 
    @time symmetrizeY!(OneLoopWorkspace,TwoLoopWorkspace,Par)

    return
end

"""Struct containing all memory used in a single ODE step """
struct MultiLoopWorkSpace{T,P <: Params{S,M} where {S,M}}
    State::StateType{T}
    Deriv::StateType{T}
    
    X::BubbleType{T} #full bubble
    XLeft::BubbleType{T} #left bubble
    XRight::BubbleType{T} #right bubble

    Y::BubbleType{T} #full bubble
    YLeft::BubbleType{T} #left bubble
    YRight::BubbleType{T} #right bubble
    Par::P
end

function Y_Workspace_Struct(Y,YTilde)
    setZero!(Y)
    setZero!(YTilde)
    return Y_Workspace_Struct(Y.x...,YTilde.x...)
end


"""
Computes a two-particle bubble in the s-Channel given two four-point functions (i.e. vertices or other bubbles).
Γ is assumed to be a vertex.
To allow the computation using just the left part of a bubble, specification of the transpose is needed i.e. Transpose(X_L) = X_R and Transpose(X) = X, where X = XL + XR is the full bubble.
"""
function compute2PartBubble!(ResultBubble,X,XTransp,Γ,getProp!,Par)
    @unpack T,N,lenIntw,np_vec,usesymmetry = Par 

    @sync begin
		for is in 1:N,it in 1:N
			Threads.@spawn begin
				BubbleProp = PropsBuffers[Threads.threadid()] # get pre-allocated thread-safe buffers
				Buffer = VertexBuffers[Threads.threadid()]
				ns = np_vec[is]
				nt = np_vec[it]
				for iu in 1:N
					nu = np_vec[iu]
					if (ns+nt+nu)%2 == 0	# skip unphysical bosonic frequency combinations
						continue
					end
					for nw in -lenIntw:lenIntw-1 # Matsubara sum
                        sprop = getProp!(BubbleProp,nw,nw+ns)
                        addYTilde!(Workspace,TwoLoopWorkspace,is,it,iu,nw,sprop,Par)
                        if(!usesymmetry || nu<=nt)
                            addY!(Workspace,TwoLoopWorkspace,is,it,iu,nw,sprop,Par,Buffer)
                        end
					end
				end
			end
		end
	end

end

"""
Computes a single-particle (i.e. self-energy) bubble.
"""
function compute1PartBubble(Result,Bubble,Prop)

end


function getTwoLoopDeriv!(Workspace::Workspace_Struct,TwoLoopWorkspace::Y_Workspace_Struct,Lam,Par,PropsBuffers,VertexBuffers)
    @unpack gamma,DVa,DVb,DVc= Workspace 
    @unpack Ya,Yb,Yc,YTa,YTb,YTc,YTd = TwoLoopWorkspace 
    @unpack T,N,Npairs,lenIntw,np_vec,usesymmetry,NUnique,OnsitePairs = Par 
    
    iG(x,nw) = iG_(gamma,x,Lam,nw,Par)
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
				Buffer = VertexBuffers[Threads.threadid()]
				ns = np_vec[is]
				nt = np_vec[it]
				for iu in 1:N
					nu = np_vec[iu]
					if (ns+nt+nu)%2 == 0	# skip unphysical bosonic frequency combinations
						continue
					end
					for nw in -lenIntw:lenIntw-1 # Matsubara sum
                        sprop = getProp!(BubbleProp,nw,nw+ns)
                        addYTilde!(Workspace,TwoLoopWorkspace,is,it,iu,nw,sprop,Par)
                        if(!usesymmetry || nu<=nt)
                            addY!(Workspace,TwoLoopWorkspace,is,it,iu,nw,sprop,Par,Buffer)
                        end
					end
				end
			end
		end
	end
end

"""Use symmetries and identities to compute the rest of bubble functions"""
function symmetrizeY!(Workspace::Workspace_Struct,TwoLoopWorkspace::Y_Workspace_Struct,Par)
    @unpack N,Npairs,usesymmetry,NUnique,OnsitePairs = Par 
    @unpack  DVa,DVb,DVc = Workspace 
    @unpack Ya,Yb,Yc,YTa,YTb,YTc,YTd = TwoLoopWorkspace 
    # use the u <--> t symmetry
    if(usesymmetry)
        Threads.@threads for it in 1:N
            for iu in it+1:N, is in 1:N, Rij in 1:Npairs
                Ya[Rij,is,it,iu] = -Ya[Rij,is,iu,it]
                Yb[Rij,is,it,iu] = -Yb[Rij,is,iu,it]
                Yc[Rij,is,it,iu] = (
                + Ya[Rij,is,it,iu]+
                - Yb[Rij,is,it,iu]+
                + Yc[Rij,is,iu,it])
            end
        end
    else
        s = 0.
        for it in 1:N, iu in it+1:N, is in 1:N, Rij in 1:Npairs
            s += abs(Ya[Rij,is,it,iu] +Ya[Rij,is,iu,it])
            s += abs( Yb[Rij,is,it,iu]  + Yb[Rij,is,iu,it])
            s += abs(Yc[Rij,is,it,iu]  -(
            + Ya[Rij,is,it,iu]+
            - Yb[Rij,is,it,iu]+
            + Yc[Rij,is,iu,it]))
        end
        println("Total Error: ",s)
    end
    #local definitions of YTilde vertices
    for iu in 1:N, it in 1:N, is in 1:N, R in OnsitePairs
        YTa[R,is,it,iu] = Ya[R,is,it,iu]
        YTb[R,is,it,iu] = Yb[R,is,it,iu]
        YTc[R,is,it,iu] = Yc[R,is,it,iu]
        YTd[R,is,it,iu] = -Yc[R,is,iu,it]
    end
    @. YTd= YTa - YTb - YTc
    Threads.@threads for iu in 1:N
        for it in 1:N, is in 1:N, Rij in 1:Npairs
            DVa[Rij,is,it,iu] += Ya[Rij,is,it,iu] - YTa[Rij,it,is,iu] + YTa[Rij,iu,is,it]
            DVb[Rij,is,it,iu] += Yb[Rij,is,it,iu] - YTc[Rij,it,is,iu] + YTc[Rij,iu,is,it]
            DVc[Rij,is,it,iu] += Yc[Rij,is,it,iu] - YTb[Rij,it,is,iu] + YTd[Rij,iu,is,it]
        end
    end
	for iu in 1:N, it in 1:N, is in 1:N, R in Par.OnsitePairs
		DVc[R,is,it,iu] = -DVb[R,it,is,iu]
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

@inline function X_(X::AbstractArray,XTransp::AbstractArray, Rj::Integer, ns::Integer,nt::Integer,nu::Integer,Rji::Integer,N::Integer)
    @assert (ns+nt+nu) %2 != 0 "$ns + $nt +  $nu = $(ns+nt+nu)"
    ns,nt,nu,swapsites = convertFreqArgs(ns,nt,nu,N)
    Vertex = ifelse(swapsites,XTransp,X)
    Rj = ifelse(swapsites,Rji,Rj)
    return @inbounds Vertex[Rj,ns+1,nt+1,nu+1]
end

@inline function XT_(XT::AbstractArray,XTTransp::AbstractArray, Rj::Integer, ns::Integer,nt::Integer,nu::Integer,Rji::Integer,N::Integer)
    @assert (ns+nt+nu) %2 != 0 "$ns + $nt +  $nu = $(ns+nt+nu)"
    Vertex = ifelse(it*iu<0,XTTransp,XT)
    ns,nt,nu,swapsites = convertFreqArgsXT(ns,nt,nu,N)
    Rj = ifelse(swapsites,Rji,Rj)
    return @inbounds Vertex[Rj,ns+1,nt+1,nu+1]
end

@inline function bufferXT_!(Cache, X::AbstractArray, XTransp::AbstractArray, ns::Integer,nt::Integer,nu::Integer,invpairs::AbstractArray,N)
    
    Vertex = ifelse(it*iu<0,XTransp,X)
    ns,nt,nu,swapsites = convertFreqArgsXT(ns,nt,nu,N)
    @assert (ns+nt+nu) %2 != 0 "$ns + $nt +  $nu = $(ns+nt+nu)"

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


struct VertexBuffer
	Va12::Vector{double}
	Vb12::Vector{double}
	Vc12::Vector{double}

	Va34::Vector{double}
	Vb34::Vector{double}
	Vc34::Vector{double}
	
	Vc21::Vector{double}
	Vc43::Vector{double}

    XTa21::Vector{double}
    XTa43::Vector{double}

    XTb21::Vector{double}
    XTb43::Vector{double}

    XTc21::Vector{double}
    XTc43::Vector{double}

    XTd21::Vector{double}
    XTd43::Vector{double}

end
VertexBuffer(Npairs) = VertexBuffer((zeros(Npairs) for _ in 1:16)...)
"""
adds to ResultBubble given the vertex from State as well as a bubble inserted on the left. XL and XR need to be provided for this.
"""
@inline function addBL!(B::BubbleType,XL::BubbleType,XR::BubbleType,State::StateType, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Props,Par::Params,Buffer)
    @unpack Va,Vb,Vc = State

    @unpack Npairs,Nsum,S,invpairs,N,np_vec = Par

    @unpack Va34,Vb34,Vc34,Vc43,XTa21,XTb21,XTc21,XTd21 = Buffer
    ns = np_vec[is]
	nt = np_vec[it]
	nu = np_vec[iu]
	wpw1,wpw2,wmw3,wmw4 = mixedFrequencies(ns,nt,nu,nwpr)

	bufferV_!(Va34, Va , ns, wmw3, wmw4, invpairs, N)
	bufferV_!(Vb34, Vb , ns, wmw3, wmw4, invpairs, N)
	bufferV_!(Vc34, Vc , ns, wmw3, wmw4, invpairs, N)
	
	bufferV_!(Vc21, Vc , ns, wpw2, wpw1, invpairs, N)
	bufferV_!(Vc43, Vc , ns, wmw4, wmw3, invpairs, N)

	bufferXT_!(XTa43, XR.Ta , XL.Ta , wmw4, ns, wmw3, invpairs, N)
	bufferXT_!(XTb43, XR.Tb , XL.Tb , wmw4, ns, wmw3, invpairs, N)
	bufferXT_!(XTc43, XR.Tc , XL.Tc , wmw4, ns, wmw3, invpairs, N)
	bufferXT_!(XTd43, XR.Td , XL.Td , wmw4, ns, wmw3, invpairs, N)

    S_ki = S.ki
	S_kj = S.kj
	S_xk = S.xk
	S_m = S.m

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
##

function addBLTilde!(B::BubbleType,XL::BubbleType,XR::BubbleType,State::StateType, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Props,Par::Params)
    @unpack Va,Vb,Vc,Xa,Xb,Xc,XTa,XTb,XTc,XTd = Workspace 
    @unpack YTa,YTb,YTc,YTd = TwoLoopWorkspace
    @unpack Npairs,invpairs,PairTypes,N,np_vec = Par
    Va_(Rij,s,t,u) = V_(Va,Rij,s,t,u,invpairs[Rij],N)
    Vb_(Rij,s,t,u) = V_(Vb,Rij,s,t,u,invpairs[Rij],N)
    Vc_(Rij,s,t,u) = V_(Vc,Rij,s,t,u,invpairs[Rij],N)

    XLa_(Rij,s,t,u) = X_(XL.a,XR.a,Rij,s,t,u,invpairs[Rij],N)
    XLb_(Rij,s,t,u) = X_(XL.b,XR.b,Rij,s,t,u,invpairs[Rij],N)
    XLc_(Rij,s,t,u) = X_(XL.c,XR.c,Rij,s,t,u,invpairs[Rij],N)

    XRTa_(Rij,s,t,u) = XT_(XR.Ta,XL.Ta,Rij,s,t,u,invpairs[Rij],N)
    XRTb_(Rij,s,t,u) = XT_(XR.Tb,XL.Tb,Rij,s,t,u,invpairs[Rij],N)
    XRTc_(Rij,s,t,u) = XT_(XR.Tc,XL.Tc,Rij,s,t,u,invpairs[Rij],N)
    XRTd_(Rij,s,t,u) = XT_(XR.Td,XL.Td,Rij,s,t,u,invpairs[Rij],N)

	ns = np_vec[is]
	nt = np_vec[it]
	nu = np_vec[iu]
	wpw1,wpw2,wmw3,wmw4 = mixedFrequencies(ns,nt,nu,nwpr)

    #Ytilde only defined for nonlocal pairs Rij >= 2
    for Rij in 2:Npairs
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