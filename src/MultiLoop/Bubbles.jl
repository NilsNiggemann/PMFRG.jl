
"""
adds to ResultBubble given the vertex as well as a bubble inserted on the left. Assumes that vertices and bubbles are given pre-computed in VertexBuffer.
"""
@inline function addBL!(B::BubbleType,XL::BubbleType,XR::BubbleType,Γ::VertexType,is::Integer, it::Integer, iu::Integer, nwpr::Integer,Par::PMFRGParams,Props,VBuffer::VertexBufferType,XBuffer::BubbleBufferType)
    @unpack N,np_vec = Par.NumericalParams
    @unpack Npairs,Nsum,siteSum,invpairs = Par.System

    fillBufferL!(VBuffer,XBuffer,XL,XR,Γ,is,it,iu,nwpr,Par)

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


            Bc_sum += (
                -Vc43[kj] * XTb21[ki] + 
                Vc34[kj] * XTd21[ki]
            )* Ptm
        end
        B.a[Rij,is,it,iu] += Ba_sum
        B.b[Rij,is,it,iu] += Bb_sum
        B.c[Rij,is,it,iu] += Bc_sum
    end
    return
end

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
"""
adds to ResultBubble given the vertex as well as a bubble inserted on the left. Assumes that vertices and bubbles are given pre-computed in VertexBuffer.
"""
@inline function addBR!(B::BubbleType,XL::BubbleType,XR::BubbleType,Γ::VertexType,is::Integer, it::Integer, iu::Integer, nwpr::Integer,Par::PMFRGParams,Props,VBuffer::VertexBufferType,XBuffer::BubbleBufferType)
    @unpack N,np_vec = Par.NumericalParams
    @unpack Npairs,Nsum,siteSum,invpairs = Par.System

    fillBufferR_new!(VBuffer,XBuffer,XL,XR,Γ,is,it,iu,nwpr,Par)

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
        Rji = invpairs[Rij]
        @turbo unroll = 1 for k_spl in 1:Nsum[Rij]
            #loop over all Nsum summation elements defined in geometry. This inner loop is responsible for most of the computational effort! 
            ki,kj,m,xk = S_ki[k_spl,Rji],S_kj[k_spl,Rji],S_m[k_spl,Rji],S_xk[k_spl,Rji]
            
            
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


            Bc_sum += (
                -Vc43[kj] * XTb21[ki] + 
                Vc34[kj] * XTd21[ki]
            )* Ptm
        end
        B.a[Rij,is,it,iu] += Ba_sum
        B.b[Rij,is,it,iu] += Bb_sum
        B.c[Rij,is,it,iu] += Bc_sum
    end
    return
end

function addBRTilde!(B::BubbleType,XL::BubbleType,XR::BubbleType,Γ::VertexType, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Par::PMFRGParams,Props)
    
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
    wpw1,wpw2,wmw3,wmw4 = wmw3,wmw4,wpw1,wpw2 #swap frequencies for right bubble
    #Btilde only defined for nonlocal pairs Rij != Rii
    for Rij in 1:Npairs
        Rij in OnsitePairs && continue
        #loop over all left hand side inequivalent pairs Rij
        Rji = invpairs[Rij] # store pair corresponding to Rji (easiest case: Rji = Rij)
        
        @unpack xi,xj = PairTypes[Rij]

        Rij, Rji = Rji, Rij #swap sites for right bubble
        xi,xj = xj,xi

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

@inline function addBL!(B::BubbleType,Γ0::BareVertexType,Γ::VertexType,is::Integer, it::Integer, iu::Integer, nwpr::Integer,Par::PMFRGParams,Props,Buffer::VertexBufferType)
    @unpack N,np_vec = Par.NumericalParams
    @unpack Npairs,Nsum,siteSum,invpairs = Par.System

    fillBufferL!(Buffer,Γ0,Γ,is,it,iu,nwpr,Par)
    
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
        Props[xj, xi]*(Va_(Rij, wmw4, ns, wmw3) + Vc_(Rij, wmw4, ns, wmw3))*Γ0.c[Rij] + 
        Props[xi, xj]*(Va_(Rji, wmw3, ns, wmw4) + Vc_(Rji, wmw3, ns, wmw4))*Γ0.c[Rji]

        B.Tc[Rij,is,it,iu] += 
        Props[xj, xi]*Vc_(Rij, wmw4, wmw3, ns)*Γ0.c[Rij] +
        Props[xi, xj]*Vc_(Rji, wmw3, wmw4, ns)*Γ0.c[Rji]
    end
end
@inline addBLTilde!(B::BubbleType,Γ0L::BareVertexType,Γ0R::BareVertexType,Γ::VertexType, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Par::PMFRGParams,Props) = addBLTilde!(B,Γ0L,Γ, is, it, iu, nwpr, Par, Props)