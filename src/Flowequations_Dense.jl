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
	Va_(Rij,s,t,u) = V_(Va,Rij,s,t,u,invpairs[Rij],N)
	Vb_(Rij,s,t,u) = V_(Vb,Rij,s,t,u,invpairs[Rij],N)

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
					jsum += (Va_(ik,0,w1pw,w1mw)+2*Vb_(ik,0,w1pw,w1mw))*iS(xk,nw)*m
				end
				Dgamma[x,iw1] += T/2 *jsum
			end
		end
	end
end


function getVertexDeriv!(Workspace::Workspace_Struct,Lam,Par,PropsBuffers,VertexBuffers)
	@unpack gamma,Dgamma,DVa,DVb,DVc,Xa,Xb,Xc,XTa,XTb,XTc,XTd = Workspace 
    @unpack T,N,Npairs,lenIntw,np_vec,usesymmetry,NUnique,OnsitePairs = Par 
	iS(x,nw) = iS_(gamma,x,Lam,nw,Par)
	iG(x,nw) = iG_(gamma,x,Lam,nw,Par)
	iSKat(x,nw) = iSKat_(gamma,Dgamma,x,Lam,nw,Par)
	function getKataninProp!(BubbleProp,nw1,nw2)
		for i in 1:NUnique, j in 1:NUnique
			BubbleProp[i,j] = iSKat(i,nw1) *iG(j,nw2)* T
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
						sprop = getKataninProp!(BubbleProp,nw,nw+ns) 
						addXTilde!(Workspace,is,it,iu,nw,sprop,Par) # add to XTilde-type bubble functions
						if(!usesymmetry || nu<=nt)
							addX!(Workspace,is,it,iu,nw,sprop,Par,Buffer)# add to X-type bubble functions
						end
					end
				end
			end
		end
	end
end

"""Use symmetries and identities to compute the rest of bubble functions"""
function symmetrizeX!(Workspace::Workspace_Struct,Par)
	@unpack DVa,DVb,DVc,Xa,Xb,Xc,XTa,XTb,XTc,XTd = Workspace 
    @unpack N,Npairs,usesymmetry,NUnique,OnsitePairs = Par 

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
	for iu in 1:N, it in 1:N, is in 1:N, R in OnsitePairs
		XTa[R,is,it,iu] = Xa[R,is,it,iu]
		XTb[R,is,it,iu] = Xb[R,is,it,iu]
		XTc[R,is,it,iu] = Xc[R,is,it,iu]
		XTd[R,is,it,iu] = -Xc[R,is,iu,it]
	end
	# SO(3) symmetry for XTd
	@. XTd= XTa - XTb - XTc
	Threads.@threads for iu in 1:N
		for it in 1:N, is in 1:N, Rij in 1:Npairs
			DVa[Rij,is,it,iu] = Xa[Rij,is,it,iu] - XTa[Rij,it,is,iu] + XTa[Rij,iu,is,it]
			DVb[Rij,is,it,iu] = Xb[Rij,is,it,iu] - XTc[Rij,it,is,iu] + XTc[Rij,iu,is,it]
			DVc[Rij,is,it,iu] = Xc[Rij,is,it,iu] - XTb[Rij,it,is,iu] + XTd[Rij,iu,is,it]
		end
	end

	for iu in 1:N, it in 1:N, is in 1:N, R in Par.OnsitePairs
		DVc[R,is,it,iu] = -DVb[R,it,is,iu]
	end
end

function mixedFrequencies(ns,nt,nu,nwpr)
	nw1=Int((ns+nt+nu-1)/2)
    nw2=Int((ns-nt-nu-1)/2)
    nw3=Int((-ns+nt-nu-1)/2)
    nw4=Int((-ns-nt+nu-1)/2)
	wpw1 = nwpr + nw1+1
    wpw2 = nwpr + nw2+1
    wmw3 = nwpr - nw3
    wmw4 = nwpr - nw4
	@assert (ns + wmw3 +wmw4)%2 != 0 "error in freq"
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
	@unpack N,Npairs,Nsum,siteSum,invpairs,np_vec = Par
	ns = np_vec[is]
	nt = np_vec[it]
	nu = np_vec[iu]
	wpw1,wpw2,wmw3,wmw4 = mixedFrequencies(ns,nt,nu,nwpr)

	bufferV_!(Va12, Va , ns, wpw1, wpw2, invpairs,N)
	bufferV_!(Vb12, Vb , ns, wpw1, wpw2, invpairs,N)
	bufferV_!(Vc12, Vc , ns, wpw1, wpw2, invpairs,N)

	bufferV_!(Va34, Va , ns, wmw3, wmw4, invpairs,N)
	bufferV_!(Vb34, Vb , ns, wmw3, wmw4, invpairs,N)
	bufferV_!(Vc34, Vc , ns, wmw3, wmw4, invpairs,N)
	
	bufferV_!(Vc21, Vc , ns, wpw2, wpw1, invpairs,N)
	bufferV_!(Vc43, Vc , ns, wmw4, wmw3, invpairs,N)
	# get fields of siteSum struct as Matrices for better use of LoopVectorization
	S_ki = siteSum.ki
	S_kj = siteSum.kj
	S_xk = siteSum.xk
	S_m = siteSum.m
	for Rij in 1:Npairs
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
	@unpack N,Npairs,invpairs,PairTypes,np_vec = Par
	Va_(Rij,s,t,u) = V_(Va,Rij,s,t,u,invpairs[Rij],N)
	Vb_(Rij,s,t,u) = V_(Vb,Rij,s,t,u,invpairs[Rij],N)
	Vc_(Rij,s,t,u) = V_(Vc,Rij,s,t,u,invpairs[Rij],N)
	ns = np_vec[is]
	nt = np_vec[it]
	nu = np_vec[iu]
	wpw1,wpw2,wmw3,wmw4 = mixedFrequencies(ns,nt,nu,nwpr)

	#Xtilde only defined for nonlocal pairs Rij >= 2
	for Rij in 2:Npairs
		#loop over all left hand side inequivalent pairs Rij
		Rji = invpairs[Rij] # store pair corresponding to Rji (easiest case: Rji = Rij)
		@unpack xi,xj = PairTypes[Rij]

		#These values are used several times so they are saved locally
		Va12 = Va_(Rji, wpw1, ns, wpw2)
		Va21 = Va_(Rij, wpw2, ns, wpw1)
		Va34 = Va_(Rji, wmw3, ns, wmw4)
		Va43 = Va_(Rij, wmw4, ns, wmw3)

		Vb12 = Vb_(Rji, wpw1, ns, wpw2)
		Vb21 = Vb_(Rij, wpw2, ns, wpw1)
		Vb34 = Vb_(Rji, wmw3, ns, wmw4)
		Vb43 = Vb_(Rij, wmw4, ns, wmw3)

		Vc12 = Vc_(Rji, wpw1, ns, wpw2)
		Vc21 = Vc_(Rij, wpw2, ns, wpw1)
		Vc34 = Vc_(Rji, wmw3, ns, wmw4)
		Vc43 = Vc_(Rij, wmw4, ns, wmw3)

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
		Vb12 = Vb_(Rji, wpw1, wpw2, ns)
		Vb21 = Vb_(Rij, wpw2, wpw1, ns)
		Vb34 = Vb_(Rji, wmw3, wmw4, ns)
		Vb43 = Vb_(Rij, wmw4, wmw3, ns)

		Vc12 = Vc_(Rji, wpw1, wpw2, ns)
		Vc21 = Vc_(Rij, wpw2, wpw1, ns)
		Vc34 = Vc_(Rji, wmw3, wmw4, ns)
		Vc43 = Vc_(Rij, wmw4, wmw3, ns)


	    XTc[Rij,is,it,iu] += (
			(+Vb21 * Vb43
			+Vc21 * Vc43
			) * Props[xi,xj]
			+(Vb12 * Vb34
			+Vc12 * Vc34
	    	)* Props[xj,xi]
		)
    end
end
##

function getChi(State, Lam::double,Par::Params,Numax)
	@unpack T,N,Npairs,lenIntw_acc,np_vec,invpairs,PairTypes,OnsitePairs = Par
	gamma = State.x[2]
	Vc = State.x[5]

	iG(x,w) = iG_(gamma,x, Lam,w,Par)
	Vc_(Rij,s,t,u) = V_(Vc,Rij,s,t,u,invpairs[Rij],N)

	Chi = zeros(Npairs,N)

	@inbounds Threads.@threads for Rij in 1:Npairs
		@unpack xi,xj = PairTypes[Rij]
		for i_nu in 1:Numax
       		n_nu = np_vec[i_nu]
		
			for nK in -lenIntw_acc:lenIntw_acc-1
				if Rij in OnsitePairs
					Chi[Rij,i_nu] += T * iG(xi,nK) * iG(xi,nK+n_nu)
				end
				for nK2 in -lenIntw_acc:lenIntw_acc-1
					npwpw2 = n_nu+nK+nK2+1
					wmw2 = nK-nK2
					#use that Vc_0 is calculated from Vb
					GGGG = iG(xi,nK)*iG(xi,nK+n_nu) * iG(xj,nK2)*iG(xj,nK2+n_nu)
					Chi[Rij,i_nu] += T^2 * GGGG *Vc_(Rij,n_nu,npwpw2,wmw2)
                end
            end
        end
    end
	return(Chi)
end

function getChi(State, Lam::double,Par::Params)
	@unpack T,N,Npairs,lenIntw_acc,np_vec,invpairs,PairTypes,OnsitePairs = Par
	gamma = State.x[2]
	Vc = State.x[5]

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
