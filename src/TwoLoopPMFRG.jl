"""includes Two-loop corrections to PMFRG. Is meant to replace FT-MFRG.jl. Do not include both in the main program
"""

include("FT-MFRG.jl")

"""Overwrite setup function"""
function setupSystem(Par::Params)
    @unpack N,Ngamma,Npairs,VDims,couplings,T,NUnique = Par
    println("TwoLoop: T= ",T)
    ##Allocate Memory:
    State = ArrayPartition(
        zeros(double,NUnique), # f_int 
        zeros(double,NUnique,Ngamma), # gamma
        zeros(double,VDims), #Va
        zeros(double,VDims), #Vb
        zeros(double,VDims) #Vc
    )

    X = CreateX(3,VDims)
    XTilde = CreateX(4,VDims)
    
    Y = CreateX(3,VDims)
    YTilde = CreateX(4,VDims)

    Vc = State.x[5]
    for is in 1:N, it in 1:N, iu in 1:N, Rj in 1:Npairs
        Vc[Rj,is,it,iu] = -couplings[Rj]
    end
    return State,(X,XTilde,Y,YTilde,Par)
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

"""Overwrite getDeriv function"""
function getDeriv!(Deriv,State,setup,Lam)
    
    X,XTilde,Y,YTilde,Par = setup #use pre-allocated X and XTilde to reduce garbage collector time
    N = Par.N
    OneLoopWorkspace = Workspace_Struct(Deriv,State,X,XTilde)
    TwoLoopWorkspace = Y_Workspace_Struct(Y,YTilde)
    @unpack DVc,DVb = OneLoopWorkspace
    getDFint!(OneLoopWorkspace,Lam,Par)
    get_Self_Energy!(OneLoopWorkspace,Lam,Par)
    getVertexDeriv!(OneLoopWorkspace,Lam,Par)
    getTwoLoopDeriv!(OneLoopWorkspace,TwoLoopWorkspace,Lam,Par)

    for iu in 1:N, it in 1:N, is in 1:N
        DVc[1,is,it,iu] = -DVb[1,it,is,iu]
    end
    return
end

function getTwoLoopDeriv!(Workspace::Workspace_Struct,TwoLoopWorkspace::Y_Workspace_Struct,Lam,Par)
	@unpack gamma,DVa,DVb,DVc= Workspace 
	@unpack Ya,Yb,Yc,YTa,YTb,YTc,YTd = TwoLoopWorkspace 
    @unpack T,N,Npairs,lenIntw,np_vec,usesymmetry,NUnique = Par 
	
	iG(x,nw) = Main.iG(gamma,x,Lam,nw,Par)

	Props = [Matrix{double}(undef,NUnique,NUnique) for _ in 1:Threads.nthreads()] 
	function getProp!(BubbleProp,nw1,nw2)
		for i in 1:NUnique, j in 1:NUnique
			BubbleProp[i,j] = iG(i,nw1) *iG(j,nw2)* T
		end
		return BubbleProp
	end
    Threads.@threads for is in 1:N
		BubbleProp = Props[Threads.threadid()]
        ns = np_vec[is]
        for nw in -lenIntw:lenIntw-1
            sprop = getProp!(BubbleProp,nw,nw+ns)
            for it in 1:N, iu in 1:N
                nt = np_vec[it]
                nu = np_vec[iu]
                if (ns+nt+nu)%2 == 0
                    continue
                end
                addYTilde!(Workspace,TwoLoopWorkspace,is,it,iu,nw,sprop,Par)
                if(!usesymmetry || nu<=nt)
                    addY!(Workspace,TwoLoopWorkspace,is,it,iu,nw,sprop,Par)
                end
            end
        end
    end
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
    end
	#local definitions of YTilde vertices
	for iu in 1:N, it in 1:N, is in 1:N
		YTa[1,is,it,iu] = Ya[1,is,it,iu]
		YTb[1,is,it,iu] = Yb[1,is,it,iu]
		YTc[1,is,it,iu] = Yc[1,is,it,iu]
		YTd[1,is,it,iu] = -Yc[1,is,iu,it]
	end
    Threads.@threads for iu in 1:N
		for it in 1:N, is in 1:N, Rij in 1:Npairs
			DVa[Rij,is,it,iu] += Ya[Rij,is,it,iu] - YTa[Rij,it,is,iu] + YTa[Rij,iu,is,it]
			DVb[Rij,is,it,iu] += Yb[Rij,is,it,iu] - YTc[Rij,it,is,iu] + YTc[Rij,iu,is,it]
			DVc[Rij,is,it,iu] += Yc[Rij,is,it,iu] - YTb[Rij,it,is,iu] + YTd[Rij,iu,is,it]
		end
	end
end

function XT_(Vertex::AbstractArray, Rj::Integer, is::Integer,it::Integer,iu::Integer,Rji::Integer)
    if is*iu<0
        return Vertex[Rji,abs(is),abs(it),abs(iu)]
    end
    return Vertex[Rj,abs(is),abs(it),abs(iu)]
end

"""
like addX! but with Y
"""
@inline function addY!(Workspace::Workspace_Struct,TwoLoopWorkspace::Y_Workspace_Struct, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Props,Par::Params)
	@unpack Va,Vb,Vc,Xa,Xb,Xc,XTa,XTb,XTc,XTd  = Workspace 
	@unpack Ya,Yb,Yc = TwoLoopWorkspace 
	@unpack Npairs,Nsum,siteSum,invpairs = Par
	
	wpw1,wpw2,wmw3,wmw4 = mixedFrequencies(is,it,iu,nwpr,Par)
    
    inverse_ki,wpw1,wpw2 = flip_tu(wpw1,wpw2) #ki always appears with w1 and w2 .. In case of XTilde, these frequencies are at first and third pos. The is freq never matters below
	inverse_kj,wmw3,wmw4 = flip_tu(wmw3,wmw4)
	@inbounds begin 
		for Rij in 1:Npairs
			#loop over all left hand side inequivalent pairs Rij
			Ya_sum = 0. #Perform summation on this temp variable before writing to State array as Base.setindex! proved to be a bottleneck!
			Yb_sum = 0.
			Yc_sum = 0.
			@fastmath @simd for k_spl in 1:Nsum[Rij]
				#loop over all Nsum summation elements defined in geometry. This inner loop is responsible for most of the computational effort! 
				@unpack ki,kj,m,xk = siteSum[k_spl,Rij]

				if (inverse_ki) 
					ki= invpairs[ki]
				end

				if (inverse_kj) 
					kj= invpairs[kj]
				end
				Ptm = Props[xk,xk]*m
				
                Ya_sum += (
                    Va[kj, is, wmw3, wmw4] * XTa[ki, wpw2, is, wpw1] + 
                    Va[ki, is, wpw1, wpw2] * XTa[kj, wmw4, is, wmw3] + 
                    2*Vb[kj, is, wmw3, wmw4] * XTc[ki, wpw2, is, wpw1] + 
                    2*Vb[ki, is, wpw1, wpw2] * XTc[kj, wmw4, is, wmw3]
				)* Ptm

				Yb_sum += (
                    Vb[kj, is, wmw3, wmw4] * XTa[ki, wpw2, is, wpw1] + 
                    Vb[ki, is, wpw1, wpw2] * XTa[kj, wmw4, is, wmw3] + 
                    Va[kj, is, wmw3, wmw4] * XTc[ki, wpw2, is, wpw1] + 
                    Vb[kj, is, wmw3, wmw4] * XTc[ki, wpw2, is, wpw1] + 
                    Va[ki, is, wpw1, wpw2] * XTc[kj, wmw4, is, wmw3] + 
                    Vb[ki, is, wpw1, wpw2] * XTc[kj, wmw4, is, wmw3]
				)* Ptm
				
				Yc_sum += (
                    -Vc[kj, is, wmw4, wmw3] * XTb[ki, wpw2, is, wpw1] - 
                    Vc[ki, is, wpw2, wpw1] * XTb[kj, wmw4, is, wmw3] + 
                    Vc[kj, is, wmw3, wmw4] * XTd[ki, wpw2, is, wpw1] + 
                    Vc[ki, is, wpw1, wpw2] * XTd[kj, wmw4, is, wmw3]
                )* Ptm
			end
			Ya[Rij,is,it,iu] += Ya_sum
			Yb[Rij,is,it,iu] += Yb_sum
			Yc[Rij,is,it,iu] += Yc_sum
        end
    end
    return
end
##

function addYTilde!(Workspace::Workspace_Struct,TwoLoopWorkspace::Y_Workspace_Struct, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Props,Par::Params)
	@unpack Va,Vb,Vc,Xa,Xb,Xc,XTa,XTb,XTc,XTd = Workspace 
    @unpack YTa,YTb,YTc,YTd = TwoLoopWorkspace
	@unpack Npairs,invpairs,PairTypes = Par
    
	Va_(Rij,s,t,u) = V_(Va,Rij,s,t,u,invpairs[Rij])
	Vb_(Rij,s,t,u) = V_(Vb,Rij,s,t,u,invpairs[Rij])
	Vc_(Rij,s,t,u) = V_(Vc,Rij,s,t,u,invpairs[Rij])

    Xa_(Rij,s,t,u) = V_(Xa,Rij,s,t,u,invpairs[Rij])
	Xb_(Rij,s,t,u) = V_(Xb,Rij,s,t,u,invpairs[Rij])
	Xc_(Rij,s,t,u) = V_(Xc,Rij,s,t,u,invpairs[Rij])

    XTa_(Rij,s,t,u) = XT_(XTa,Rij,s,t,u,invpairs[Rij])
	XTb_(Rij,s,t,u) = XT_(XTb,Rij,s,t,u,invpairs[Rij])
	XTc_(Rij,s,t,u) = XT_(XTc,Rij,s,t,u,invpairs[Rij])
	XTd_(Rij,s,t,u) = XT_(XTd,Rij,s,t,u,invpairs[Rij])

	wpw1,wpw2,wmw3,wmw4 = mixedFrequencies(is,it,iu,nwpr,Par)

	#Ytilde only defined for nonlocal pairs Rij >= 2
	for Rij in 2:Npairs
		#loop over all left hand side inequivalent pairs Rij
		Rji = invpairs[Rij] # store pair corresponding to Rji (easiest case: Rji = Rij)
		@unpack xi,xj = PairTypes[Rij]

	    YTa[Rij,is,it,iu] += 
        Props[xj, xi]*(
        Va_(Rij, wpw2, is, wpw1)*Xa_(Rij, wmw4, is, wmw3) + 
        Va_(Rij, wmw4, is, wmw3)*Xa_(Rij, wpw2, is, wpw1) + 
        2*Vc_(Rij, wpw2, is, wpw1)*Xc_(Rij, wmw4, is, wmw3) + 
        2*Vc_(Rij, wmw4, is, wmw3)*Xc_(Rij, wpw2, is, wpw1)) + 
        Props[xi, xj]*(
        Va_(Rji, wpw1, is, wpw2)*XTa_(Rji, wmw4, wmw3, is) + 
        Va_(Rji, wmw3, is, wmw4)*XTa_(Rji, wpw2, wpw1, is) + 
        2*Vc_(Rji, wpw1, is, wpw2)*XTd_(Rji, wmw4, wmw3, is) + 
        2*Vc_(Rji, wmw3, is, wmw4)*XTd_(Rji, wpw2, wpw1, is))
		
	    YTb[Rij,is,it,iu] += 
        Props[xj, xi]*(
        Va_(Rij, wpw2, is, wpw1)*Xc_(Rij, wmw4, is, wmw3) + 
        Vc_(Rij, wpw2, is, wpw1)*(
        Xa_(Rij, wmw4, is, wmw3) + Xc_(Rij, wmw4, is, wmw3)) + 
        Va_(Rij, wmw4, is, wmw3)*Xc_(Rij, wpw2, is, wpw1) + 
        Vc_(Rij, wmw4, is, wmw3)*(
        Xa_(Rij, wpw2, is, wpw1) + Xc_(Rij, wpw2, is, wpw1))) + 
        Props[xi, xj]*(
        Va_(Rji, wpw1, is, wpw2)*XTd_(Rji, wmw4, wmw3, is) + 
        Vc_(Rji, wpw1, is, wpw2)*(
        XTa_(Rji, wmw4, wmw3, is) + XTd_(Rji, wmw4, wmw3, is)) + 
        Va_(Rji, wmw3, is, wmw4)*XTd_(Rji, wpw2, wpw1, is) + 
        Vc_(Rji, wmw3, is, wmw4)*(
        XTa_(Rji, wpw2, wpw1, is) + XTd_(Rji, wpw2, wpw1, is)))

	    YTc[Rij,is,it,iu] += 
        Props[xj, xi]*(
        Vb_(Rij, wpw2, is, wpw1)*Xb_(Rij, wmw4, is, wmw3) + 
        Vb_(Rij, wmw4, is, wmw3)*Xb_(Rij, wpw2, is, wpw1) + 
        Vc_(Rij, wpw2, wpw1, is)*Xc_(Rij, wmw4, wmw3, is) + 
        Vc_(Rij, wmw4, wmw3, is)*Xc_(Rij, wpw2, wpw1, is)) + 
        Props[xi, xj]*(
        -Vc_(Rji, wpw1, wpw2, is)*XTb_(Rji, wmw4, wmw3, is) - 
        Vc_(Rji, wmw3, wmw4, is)*XTb_(Rji, wpw2, wpw1, is) + 
        Vb_(Rji, wpw1, is, wpw2)*XTc_(Rji, wmw4, wmw3, is) + 
        Vb_(Rji, wmw3, is, wmw4)*XTc_(Rji, wpw2, wpw1, is))

	    YTd[Rij,is,it,iu] += 
        -Props[xj, xi]*(
        Vc_(Rij, wpw2, wpw1, is)*Xb_(Rij, wmw4, is, wmw3) + 
        Vc_(Rij, wmw4, wmw3, is)*Xb_(Rij, wpw2, is, wpw1) + 
        Vb_(Rij, wpw2, is, wpw1)*Xc_(Rij, wmw4, wmw3, is) + 
        Vb_(Rij, wmw4, is, wmw3)*Xc_(Rij, wpw2, wpw1, is)) + 
        Props[xi, xj]*(
        Vb_(Rji, wpw1, is, wpw2)*XTb_(Rji, wmw4, wmw3, is) + 
        Vb_(Rji, wmw3, is, wmw4)*XTb_(Rji, wpw2, wpw1, is) - 
        Vc_(Rji, wpw1, wpw2, is)*XTc_(Rji, wmw4, wmw3, is) - 
        Vc_(Rji, wmw3, wmw4, is)*XTc_(Rji, wpw2, wpw1, is))
    end
end