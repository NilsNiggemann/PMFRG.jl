function getDFint!(Workspace::Workspace_Struct,Lam::double,Par::Params)
    @unpack Df_int,gamma = Workspace 
    @unpack T,lenIntw_acc,NUnique = Par 
	
	γ(x,nw) = gamma_(gamma,x,nw,Par)
	iG(x,nw) = iG_(gamma,x,Lam,nw,Par)
	iS(x,nw) = iS_(gamma,x,Lam,nw,Par)

	Theta(Lam,w) = w^2/(w^2+Lam^2)
	
	for x in 1:NUnique
		sumres = 0.
		for nw in -lenIntw_acc:lenIntw_acc-1
			w = get_w(nw,T)
			sumres += iS(x,nw)/iG(x,nw)*Theta(Lam,w) *γ(x,nw)/w
		end
		Df_int[x] = -3/2*T*sumres
	end
end


function get_Self_Energy!(Workspace::Workspace_Struct,Lam::double,Par::Params)
    @unpack Dgamma,gamma,Va,Vb = Workspace 
    @unpack T,N,Ngamma,lenIntw_acc,np_vec_gamma,siteSum,invpairs,Nsum,OnsitePairs = Par 
	
	iS(x,nw) = iS_(gamma,x,Lam,nw,Par)
	Va_(Rij,s,t,u) = V_(Va,Rij,s,t,u,invpairs[Rij])
	Vb_(Rij,s,t,u) = V_(Vb,Rij,s,t,u,invpairs[Rij])

	Threads.@threads for iw1 in 1:Ngamma
		nw1 = np_vec_gamma[iw1]
    	for (x,Rx) in enumerate(OnsitePairs)
			for nw in -lenIntw_acc: lenIntw_acc-1
				jsum = 0.
				w1pw = get_sign_iw(nw1+nw+1,N) #w1 + w: Adding two Matsubara frequencies gives a +1
				w1mw = get_sign_iw(nw1-nw,N)
				for k_spl in 1:Nsum[Rx]
					@unpack m,ki,xk = siteSum[k_spl,Rx]
					ik = invpairs[ki] # pair ik is inversed relative to pre-generated site summation in X (ki)!
					jsum += (Va_(ik,1,w1pw,w1mw)+2*Vb_(ik,1,w1pw,w1mw))*iS(xk,nw)*m
				end
				Dgamma[x,iw1] += T/2 *jsum
			end
		end
	end
end

function getVertexDeriv!(Workspace::Workspace_Struct,Lam,Par)
	@unpack gamma,Dgamma,DVa,DVb,DVc,Xa,Xb,Xc,XTa,XTb,XTc,XTd = Workspace 
    @unpack T,N,Npairs,lenIntw,np_vec,usesymmetry,NUnique = Par 
	
	iS(x,nw) = iS_(gamma,x,Lam,nw,Par)
	iG(x,nw) = iG_(gamma,x,Lam,nw,Par)
	iSKat(x,nw) = iSKat_(gamma,Dgamma,x,Lam,nw,Par)

	Props = [Matrix{double}(undef,NUnique,NUnique) for _ in 1:Threads.nthreads()] 
	Buffers = [VertexBuffer(Par.Npairs) for _ in 1:Threads.nthreads()] 
	function getKataninProp!(BubbleProp,nw1,nw2)
		for i in 1:NUnique, j in 1:NUnique
			BubbleProp[i,j] = iSKat(i,nw1) *iG(j,nw2)* T
		end
		return BubbleProp
	end
    Threads.@threads for is in 1:N
		BubbleProp = Props[Threads.threadid()]
		Buffer = Buffers[Threads.threadid()]
        ns = np_vec[is]
        for nw in -lenIntw:lenIntw-1
            sprop = getKataninProp!(BubbleProp,nw,nw+ns)
            for it in 1:N, iu in 1:N
                nt = np_vec[it]
                nu = np_vec[iu]
                if (ns+nt+nu)%2 == 0
                    continue
                end
                addXTilde!(Workspace,is,it,iu,nw,sprop,Par)
                if(!usesymmetry || nu<=nt)
                    addX!(Workspace,is,it,iu,nw,sprop,Par,Buffer)
                end
            end
        end
    end
    # use the u <--> t symmetry
    if(usesymmetry)
		Threads.@threads for it in 1:N
			for iu in it+1:N, is in 1:N, Rij in 1:Npairs
				Xa[Rij,is,it,iu] = -Xa[Rij,is,iu,it]
				Xb[Rij,is,it,iu] = -Xb[Rij,is,iu,it]
				Xc[Rij,is,it,iu] = (
				+ Xa[Rij,is,it,iu]+
				- Xb[Rij,is,it,iu]+
				+ Xc[Rij,is,iu,it])
			end
		end
    end
	#local definitions of XTilde vertices
	for iu in 1:N, it in 1:N, is in 1:N
		XTa[1,is,it,iu] = Xa[1,is,it,iu]
		XTb[1,is,it,iu] = Xb[1,is,it,iu]
		XTc[1,is,it,iu] = Xc[1,is,it,iu]
		XTd[1,is,it,iu] = -Xc[1,is,iu,it]
	end
	Threads.@threads for iu in 1:N
		for it in 1:N, is in 1:N, Rij in 1:Npairs
			DVa[Rij,is,it,iu] = Xa[Rij,is,it,iu] - XTa[Rij,it,is,iu] + XTa[Rij,iu,is,it]
			DVb[Rij,is,it,iu] = Xb[Rij,is,it,iu] - XTc[Rij,it,is,iu] + XTc[Rij,iu,is,it]
			DVc[Rij,is,it,iu] = Xc[Rij,is,it,iu] - XTb[Rij,it,is,iu] + XTd[Rij,iu,is,it]
		end
	end

end

function mixedFrequencies(is,it,iu,nwpr,Par)
	@unpack np_vec,N = Par
    ns = np_vec[is]
	nt = np_vec[it]
	nu = np_vec[iu]

	nw1=div(ns+nt+nu+1,2)
    nw2=div(ns-nt-nu+1,2)
    nw3=div(-ns+nt-nu+1,2)
    nw4=div(-ns-nt+nu+1,2)
	wpw1 = get_sign_iw(nwpr + nw1+1,N)
    wpw2 = get_sign_iw(nwpr + nw2+1,N)
    wmw3 = get_sign_iw(nwpr - nw3,N)
    wmw4 = get_sign_iw(nwpr - nw4,N)
	return wpw1,wpw2,wmw3,wmw4
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
end
VertexBuffer(Npairs) = VertexBuffer((zeros(Npairs) for _ in 1:8)...)
"""
adds part of X functions in Matsubara sum at nwpr containing the site summation for a set of s t and u frequencies. This is the most numerically demanding part!
"""
@inline function addX!(Workspace::Workspace_Struct, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Props,Par::Params,Buffer)
	@unpack Va,Vb,Vc,Xa,Xb,Xc = Workspace 
	@unpack Va12,Vb12,Vc12,Va34,Vb34,Vc34,Vc21,Vc43 = Buffer 
	@unpack Npairs,Nsum,S,invpairs = Par

	wpw1,wpw2,wmw3,wmw4 = mixedFrequencies(is,it,iu,nwpr,Par)

	bufferV_!(Va12, Va , is, wpw1, wpw2, invpairs)
	bufferV_!(Vb12, Vb , is, wpw1, wpw2, invpairs)
	bufferV_!(Vc12, Vc , is, wpw1, wpw2, invpairs)

	bufferV_!(Va34, Va , is, wmw3, wmw4, invpairs)
	bufferV_!(Vb34, Vb , is, wmw3, wmw4, invpairs)
	bufferV_!(Vc34, Vc , is, wmw3, wmw4, invpairs)
	
	bufferV_!(Vc21, Vc , is, wpw2, wpw1, invpairs)
	bufferV_!(Vc43, Vc , is, wmw4, wmw3, invpairs)
	# get fields of siteSum struct as Matrices for better use of LoopVectorization
	S_ki = S.ki
	S_kj = S.kj
	S_xk = S.xk
	S_m = S.m
	@inbounds for Rij in 1:Npairs
		#loop over all left hand side inequivalent pairs Rij
		Xa_sum = 0. #Perform summation on this temp variable before writing to State array as Base.setindex! proved to be a bottleneck!
		Xb_sum = 0.
		Xc_sum = 0.
		@turbo unroll = 1 for k_spl in 1:Nsum[Rij]
			#loop over all Nsum summation elements defined in geometry. This inner loop is responsible for most of the computational effort! 
			# @unpack ki,kj,m,xk = siteSum[k_spl,Rij]
			ki,kj,m,xk = S_ki[k_spl,Rij],S_kj[k_spl,Rij],S_m[k_spl,Rij],S_xk[k_spl,Rij]
			Ptm = Props[xk,xk]*m

			Xa_sum += (
				+Va12[ki] * Va34[kj] 
				+Vb12[ki] * Vb34[kj] * 2
			)* Ptm

			Xb_sum += (
				+Va12[ki] * Vb34[kj]
				+Vb12[ki] * Va34[kj]
				+Vb12[ki] * Vb34[kj]
			)* Ptm
			
			Xc_sum += (
				+Vc12[ki] * Vc34[kj]
				+Vc21[ki] * Vc43[kj]
			)* Ptm
		end
		Xa[Rij,is,it,iu] += Xa_sum
		Xb[Rij,is,it,iu] += Xb_sum
		Xc[Rij,is,it,iu] += Xc_sum
    end
    return
end
##
function addXTilde!(Workspace::Workspace_Struct, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Props,Par::Params)
	@unpack Va,Vb,Vc,XTa,XTb,XTc,XTd = Workspace 
	@unpack Npairs,invpairs,PairTypes = Par
    
	Va_(Rij,s,t,u) = V_(Va,Rij,s,t,u,invpairs[Rij])
	Vb_(Rij,s,t,u) = V_(Vb,Rij,s,t,u,invpairs[Rij])
	Vc_(Rij,s,t,u) = V_(Vc,Rij,s,t,u,invpairs[Rij])

	wpw1,wpw2,wmw3,wmw4 = mixedFrequencies(is,it,iu,nwpr,Par)

	#Xtilde only defined for nonlocal pairs Rij >= 2
	for Rij in 2:Npairs
		#loop over all left hand side inequivalent pairs Rij
		Rji = invpairs[Rij] # store pair corresponding to Rji (easiest case: Rji = Rij)
		@unpack xi,xj = PairTypes[Rij]

		#These values are used several times so they are saved locally
		Va12 = Va_(Rji, wpw1, is, wpw2)
		Va21 = Va_(Rij, wpw2, is, wpw1)
		Va34 = Va_(Rji, wmw3, is, wmw4)
		Va43 = Va_(Rij, wmw4, is, wmw3)

		Vb12 = Vb_(Rji, wpw1, is, wpw2)
		Vb21 = Vb_(Rij, wpw2, is, wpw1)
		Vb34 = Vb_(Rji, wmw3, is, wmw4)
		Vb43 = Vb_(Rij, wmw4, is, wmw3)

		Vc12 = Vc_(Rji, wpw1, is, wpw2)
		Vc21 = Vc_(Rij, wpw2, is, wpw1)
		Vc34 = Vc_(Rji, wmw3, is, wmw4)
		Vc43 = Vc_(Rij, wmw4, is, wmw3)

	    XTa[Rij,is,it,iu] += (
			(+Va21 * Va43
			+2*Vc21 * Vc43) * Props[xi,xj]
			+(Va12 * Va34
			+2*Vc12 * Vc34)* Props[xj,xi]
		)
		
	    XTb[Rij,is,it,iu] += (
			(+Va21 * Vc43
			+Vc21 * Vc43
			+Vc21 * Va43) * Props[xi,xj]

			+(Va12 * Vc34
			+Vc12 * Vc34
			+Vc12 * Va34)* Props[xj,xi]
		)
		Vb12 = Vb_(Rji, wpw1, wpw2, is)
		Vb21 = Vb_(Rij, wpw2, wpw1, is)
		Vb34 = Vb_(Rji, wmw3, wmw4, is)
		Vb43 = Vb_(Rij, wmw4, wmw3, is)

		Vc12 = Vc_(Rji, wpw1, wpw2, is)
		Vc21 = Vc_(Rij, wpw2, wpw1, is)
		Vc34 = Vc_(Rji, wmw3, wmw4, is)
		Vc43 = Vc_(Rij, wmw4, wmw3, is)


	    XTc[Rij,is,it,iu] += (
			(+Vb21 * Vb43
			+Vc21 * Vc43
			) * Props[xi,xj]
			+(Vb12 * Vb34
			+Vc12 * Vc34
	    	)* Props[xj,xi]
		)

	    XTd[Rij,is,it,iu] += (
			(+Vb21 * Vc43
			+Vc21 * Vb43
			) * Props[xi,xj]
			+(Vb12 * Vc34
			+Vc12 * Vb34
	    	)* Props[xj,xi]
		)
    end
end
##

function getChi(State, Lam::double,Par::Params,Numax = 1)
	@unpack T,N,Npairs,lenIntw_acc,np_vec,invpairs,PairTypes,OnsitePairs = Par
	gamma = State.x[2]
	Vc = State.x[5]

	iG(x,w) = iG_(gamma,x, Lam,w,Par)
	Vc_(Rij,s,t,u) = V_(Vc,Rij,s,t,u,invpairs[Rij])

	Chi = zeros(Npairs,N)
	Threads.@threads for Rij in 1:Npairs
		@unpack xi,xj = PairTypes[Rij]
		for i_nu in 1:Numax
       		n_nu = np_vec[i_nu]
		
			for nK in -lenIntw_acc:lenIntw_acc-1
				if Rij in OnsitePairs
					Chi[Rij,i_nu] += T * iG(xi,nK) * iG(xi,nK+n_nu)
				end
				for nK2 in -lenIntw_acc:lenIntw_acc-1
					npwpw2 = get_sign_iw(n_nu+nK+nK2+1,N)
					wmw2 = get_sign_iw(nK-nK2,N)
					#use that Vc_0 is calculated from Vb
					GGGG = iG(xi,nK)*iG(xi,nK+n_nu) * iG(xj,nK2)*iG(xj,nK2+n_nu)
					Chi[Rij,i_nu] += T^2 * GGGG *Vc_(Rij,i_nu,npwpw2,wmw2)
                end
            end
        end
    end
	return(Chi)
end
