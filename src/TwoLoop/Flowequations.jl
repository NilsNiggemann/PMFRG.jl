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

struct Y_Workspace_Struct
    Ya::Array{double,4}
    Yb::Array{double,4}
    Yc::Array{double,4}
    YTa::Array{double,4}
    YTb::Array{double,4}
    YTc::Array{double,4}
    YTd::Array{double,4}
end

function Y_Workspace_Struct(Y,YTilde)
    setZero!(Y)
    setZero!(YTilde)
    return Y_Workspace_Struct(Y.x...,YTilde.x...)
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

function XT_(Vertex::AbstractArray, Rj::Integer, ns::Integer,nt::Integer,nu::Integer,Rji::Integer,N::Integer)
    @assert (ns+nt+nu) %2 != 0 "$ns + $nt +  $nu = $(ns+nt+nu)"
    ns,nt,nu,swapsites = convertFreqArgsXT(ns,nt,nu,N)
    Rj = ifelse(swapsites,Rji,Rj)
    return @inbounds Vertex[Rj,ns+1,nt+1,nu+1]
end

@inline function bufferXT_!(Cache, Vertex::AbstractArray, ns::Integer,nt::Integer,nu::Integer,invpairs::AbstractArray,N)
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
like addX! but with Y
"""
@inline function addY!(Workspace::Workspace_Struct,TwoLoopWorkspace::Y_Workspace_Struct, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Props,Par::Params,Buffer)
    @unpack Va,Vb,Vc,XTa,XTb,XTc,XTd  = Workspace 
    @unpack Ya,Yb,Yc = TwoLoopWorkspace 
    @unpack Npairs,Nsum,S,invpairs,N,np_vec = Par
    @unpack Va12,Vb12,Vc12,Va34,Vb34,Vc34,Vc21,Vc43,XTa21,XTa43,XTb21,XTb43,XTc21,XTc43,XTd21,XTd43 = Buffer
    ns = np_vec[is]
	nt = np_vec[it]
	nu = np_vec[iu]
	wpw1,wpw2,wmw3,wmw4 = mixedFrequencies(ns,nt,nu,nwpr)

    bufferV_!(Va12, Va , ns, wpw1, wpw2, invpairs, N)
	bufferV_!(Vb12, Vb , ns, wpw1, wpw2, invpairs, N)
	bufferV_!(Vc12, Vc , ns, wpw1, wpw2, invpairs, N)

	bufferV_!(Va34, Va , ns, wmw3, wmw4, invpairs, N)
	bufferV_!(Vb34, Vb , ns, wmw3, wmw4, invpairs, N)
	bufferV_!(Vc34, Vc , ns, wmw3, wmw4, invpairs, N)
	
	bufferV_!(Vc21, Vc , ns, wpw2, wpw1, invpairs, N)
	bufferV_!(Vc43, Vc , ns, wmw4, wmw3, invpairs, N)

    bufferXT_!(XTa21, XTa , wpw2, ns, wpw1, invpairs, N)
	bufferXT_!(XTa43, XTa , wmw4, ns, wmw3, invpairs, N)

    bufferXT_!(XTb21, XTb , wpw2, ns, wpw1, invpairs, N)
	bufferXT_!(XTb43, XTb , wmw4, ns, wmw3, invpairs, N)

    bufferXT_!(XTc21, XTc , wpw2, ns, wpw1, invpairs, N)
	bufferXT_!(XTc43, XTc , wmw4, ns, wmw3, invpairs, N)

    bufferXT_!(XTd21, XTd , wpw2, ns, wpw1, invpairs, N)
	bufferXT_!(XTd43, XTd , wmw4, ns, wmw3, invpairs, N)

    S_ki = S.ki
	S_kj = S.kj
	S_xk = S.xk
	S_m = S.m

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
        Ya[Rij,is,it,iu] += Ya_sum
        Yb[Rij,is,it,iu] += Yb_sum
        Yc[Rij,is,it,iu] += Yc_sum
    end
    return
end
##

function addYTilde!(Workspace::Workspace_Struct,TwoLoopWorkspace::Y_Workspace_Struct, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Props,Par::Params)
    @unpack Va,Vb,Vc,Xa,Xb,Xc,XTa,XTb,XTc,XTd = Workspace 
    @unpack YTa,YTb,YTc,YTd = TwoLoopWorkspace
    @unpack Npairs,invpairs,PairTypes,N,np_vec = Par
    Va_(Rij,s,t,u) = V_(Va,Rij,s,t,u,invpairs[Rij],N)
    Vb_(Rij,s,t,u) = V_(Vb,Rij,s,t,u,invpairs[Rij],N)
    Vc_(Rij,s,t,u) = V_(Vc,Rij,s,t,u,invpairs[Rij],N)

    Xa_(Rij,s,t,u) = V_(Xa,Rij,s,t,u,invpairs[Rij],N)
    Xb_(Rij,s,t,u) = V_(Xb,Rij,s,t,u,invpairs[Rij],N)
    Xc_(Rij,s,t,u) = V_(Xc,Rij,s,t,u,invpairs[Rij],N)

    XTa_(Rij,s,t,u) = XT_(XTa,Rij,s,t,u,invpairs[Rij],N)
    XTb_(Rij,s,t,u) = XT_(XTb,Rij,s,t,u,invpairs[Rij],N)
    XTc_(Rij,s,t,u) = XT_(XTc,Rij,s,t,u,invpairs[Rij],N)
    XTd_(Rij,s,t,u) = XT_(XTd,Rij,s,t,u,invpairs[Rij],N)

	ns = np_vec[is]
	nt = np_vec[it]
	nu = np_vec[iu]
	wpw1,wpw2,wmw3,wmw4 = mixedFrequencies(ns,nt,nu,nwpr)

    #Ytilde only defined for nonlocal pairs Rij >= 2
    for Rij in 2:Npairs
        #loop over all left hand side inequivalent pairs Rij
        Rji = invpairs[Rij] # store pair corresponding to Rji (easiest case: Rji = Rij)
        @unpack xi,xj = PairTypes[Rij]

        YTa[Rij,is,it,iu] += 
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
        
        YTb[Rij,is,it,iu] += 
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

        YTc[Rij,is,it,iu] += 
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