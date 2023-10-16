function getDeriv!(Deriv, State, setup::Tuple{BubbleType,Ty,OneLoopParams}, T) where {Ty}
    (X, Buffs, Par) = setup #use pre-allocated X and XTilde to reduce garbage collector time
    Workspace = OneLoopWorkspace(Deriv, State, X, Buffs, Par)

    getDFint!(Workspace, T)
    get_Self_Energy!(Workspace, T)

    getXBubble!(Workspace,T)
    symmetrizeBubble!(Workspace.X, Par)

    addToVertexFromBubble!(Workspace.Deriv.Γ, Workspace.X)
    symmetrizeVertex!(Workspace.Deriv.Γ, Par)
    # flush(stdout)
    return
end
function getXBubble!(Workspace,T)
    getX!(Workspace, T)
    getXTilde!(Workspace,T)
end

function getDFint!(Workspace::PMFRGWorkspace, T::Real)
    (; State, Deriv, Par) = Workspace
    (; lenIntw_acc) = Par.NumericalParams
    NUnique = Par.System.NUnique

    @inline γ(x, nw) = gamma_(State.γ, x, nw)
    @inline iG(x, nw) = iG_(State.γ, x, T, nw)
    @inline iS(x, nw) = iS_(State.γ, x, T, nw)

    for x = 1:NUnique
        sumres = 0.0
        for nw = -lenIntw_acc:lenIntw_acc-1
            sumres += iS(x, nw) / iG(x, nw) * iG0(nw, T) * γ(x, nw)
        end
        Deriv.f_int[x] = -3 / 2 * sumres
    end
end

"""
Computes a single-particle (i.e. self-energy) bubble. Allows specification of function type, i.e. what vertices are used since this is different if a bubble function is inserted as opposed to a vertex.
"""
function addTo1PartBubble!(Dgamma::AbstractArray, XT1_::Function, XT2_::Function, Prop, Par)

    (; N, Ngamma, lenIntw_acc, np_vec_gamma) = Par.NumericalParams
    (; siteSum, invpairs, Nsum, OnsitePairs) = Par.System

    Threads.@threads for iw1 = 1:Ngamma
        nw1 = np_vec_gamma[iw1]
        for (x, Rx) in enumerate(OnsitePairs)
            for nw = -lenIntw_acc:lenIntw_acc-1
                jsum = 0.0
                wpw1 = nw1 + nw + 1 #w + w1: Adding two fermionic Matsubara frequencies gives a +1 for the bosonic index
                wmw1 = nw - nw1
                for k_spl = 1:Nsum[Rx]
                    (; m, ki, xk) = siteSum[k_spl, Rx]
                    jsum +=
                        (XT1_(ki, wpw1, 0, wmw1) + 2 * XT2_(ki, wpw1, 0, wmw1)) *
                        Prop(xk, nw) *
                        m
                end
                Dgamma[x, iw1] += -jsum #For the self-energy derivative, the factor of 1/2 must be in the propagator
            end
        end
    end
    return Dgamma
end

"""
Computes a single-particle (i.e. self-energy) bubble. Can only be used if B is a bubble function
"""
function addTo1PartBubble!(Dgamma::AbstractArray, X::BubbleType, Prop, Par)
    invpairs = Par.System.invpairs
    @inline XTa_(Rij, s, t, u) =
        XT_(X.Ta, X.Ta, Rij, s, t, u, invpairs[Rij], Par.NumericalParams.N)
    @inline XTc_(Rij, s, t, u) =
        XT_(X.Tc, X.Tc, Rij, s, t, u, invpairs[Rij], Par.NumericalParams.N)
    addTo1PartBubble!(Dgamma, XTa_, XTc_, Prop, Par)
end

"""
Computes a single-particle (i.e. self-energy) bubble. Can only be used if argument is a vertex
"""
function addTo1PartBubble!(Dgamma::AbstractArray, Γ::VertexType, Prop, Par)
    invpairs = Par.System.invpairs
    # @warn "addTo1PartBubble! for vertices is not tested yet!"
    @inline Γa_(Rij, s, t, u) = V_(Γ.a, Rij, t, u, s, invpairs[Rij], Par.NumericalParams.N) # Tilde-type can be obtained by permutation of vertices
    @inline Γb_(Rij, s, t, u) = V_(Γ.b, Rij, t, u, s, invpairs[Rij], Par.NumericalParams.N) # cTilde corresponds to b type vertex!
    addTo1PartBubble!(Dgamma, Γa_, Γb_, Prop, Par)
end

function compute1PartBubble!(Dgamma::AbstractArray, ΓorX, Prop, Par)
    setZero!(Dgamma)
    addTo1PartBubble!(Dgamma, ΓorX, Prop, Par)
end
function get_Self_Energy!(Workspace::PMFRGWorkspace, T)
    Par = Workspace.Par
    @inline iS(x, nw) = iS_(Workspace.State.γ, x, T, nw) / 2
    compute1PartBubble!(Workspace.Deriv.γ, Workspace.State.Γ, iS, Par)
end

function getX!(Workspace::PMFRGWorkspace, T)
    Par = Workspace.Par
    (; N, np_vec) = Par.NumericalParams
    PropsBuffers = Workspace.Buffer.Props
    VertexBuffers = Workspace.Buffer.Vertex
    for is = 1:N
        ns = np_vec[is]
        for Rij = 1:Par.System.Npairs
            bufferPropagator!(PropsBuffers,Rij, ns, Workspace,T)
            bufferVertices!(VertexBuffers, Rij, is, Workspace)
            @sync begin
                for it = 1:N, iu = 1:N
                    nt = np_vec[it]
                    nu = np_vec[iu]
                    (ns + nt + nu) % 2 == 0 && continue # skip unphysical bosonic frequency combinations
                    nu >nt && Par.Options.usesymmetry && continue
                    Threads.@spawn addX!(Workspace, Rij, is, it, iu, PropsBuffers, VertexBuffers)# add to X-type bubble functions
                end
            end
        end
    end
end

function bufferPropagator!(PropBuffer,Rij,ns,Workspace,T)
    Par = Workspace.Par
    (; Nsum, siteSum) = Par.System
    setZero!(PropBuffer)
    iG(x, nw) = iG_(Workspace.State.γ, x, T, nw)
    iSKat(x, nw) = iSKat_(Workspace.State.γ, Workspace.Deriv.γ, x, T, nw)

    lenIntw = Par.NumericalParams.lenIntw

    for k_spl = 1:Nsum[Rij]
        m = siteSum.m[k_spl,Rij]
        xk = siteSum.xk[k_spl,Rij]
        for (i,nw) in enumerate(-lenIntw:lenIntw-1)
            PropBuffer[k_spl,i] = iSKat(xk, nw) * iG(xk, nw+ns) * m
        end
    end
    return PropBuffer
end

Base.@propagate_inbounds function bufferV!(V_ki_kj,Vertex,Rij,is,Par)
    (; Nsum, siteSum, invpairs) = Par.System
    setZero!(V_ki_kj)
    Base.require_one_based_indexing(V_ki_kj,Vertex)
    @boundscheck size(V_ki_kj) === (Nsum[Rij],2,2,Par.NumericalParams.N,Par.NumericalParams.N)
    @boundscheck size(V_ki_kj) === (Nsum[Rij],2,2,Par.NumericalParams.N,Par.NumericalParams.N)
    @inbounds @simd for iu in axes(Vertex,4)
        for it in axes(Vertex,3)
            for k_spl = 1:Nsum[Rij]
                ki = siteSum.ki[k_spl,Rij]
                kj = siteSum.kj[k_spl,Rij]
                ik = invpairs[ki]
                jk = invpairs[kj]
            #   V[k,side,swapsites,nt,nu]
                V_ki_kj[k_spl,1,1,it,iu] = Vertex[ki,is, it, iu] # first index is the site index, second is for either ki or kj, third is for invpairs
                V_ki_kj[k_spl,1,2,it,iu] = Vertex[ik,is, it, iu]

                V_ki_kj[k_spl,2,1,it,iu] = Vertex[kj,is, it, iu]
                V_ki_kj[k_spl,2,2,it,iu] = Vertex[jk,is, it, iu]
            end
        end
    end

    return V_ki_kj
end

function getBufferView(V_ki_kj,side,ns,nt,nu,N,Nsum)
    ns, nt, nu, swapsites = convertFreqArgs(ns, nt, nu, N)
    # @assert (ns+nt+nu) %2 != 0 "$ns + $nt +  $nu = $(ns+nt+nu)"
    it = nt + 1
    iu = nu + 1

    if swapsites
        return @view V_ki_kj[begin:Nsum,side,2,it,iu]
    else
        return @view V_ki_kj[begin:Nsum,side,1,it,iu]
    end
end

getBufferView_ki(Vki_kj,ns,nt,nu,N,Nsum) = getBufferView(Vki_kj,1,ns,nt,nu,N,Nsum)
getBufferView_kj(Vki_kj,ns,nt,nu,N,Nsum) = getBufferView(Vki_kj,2,ns,nt,nu,N,Nsum)

function bufferVertices!(VertexBuffers,Rij,is,Workspace)
    (;Par,State) = Workspace
    (;Va_ki_kj, Vb_ki_kj, Vc_ki_kj) = VertexBuffers

    bufferV!(Va_ki_kj,State.Γ.a,Rij,is,Par)
    bufferV!(Vb_ki_kj,State.Γ.b,Rij,is,Par)
    bufferV!(Vc_ki_kj,State.Γ.c,Rij,is,Par)
end



function getXTilde!(Workspace::PMFRGWorkspace, T)
    Par = Workspace.Par
    (; N, lenIntw, np_vec) = Par.NumericalParams
    @sync begin
        for is = 1:N
            ns = np_vec[is]
            for it = 1:N
                nt = np_vec[it]
                Threads.@spawn begin
                    # Workspace.X.a .= Buffer.Va12[begin]
                    for iu = 1:N
                        nu = np_vec[iu]
                        (ns + nt + nu) % 2 == 0 && continue
                        for nw = -lenIntw:lenIntw-1 # Matsubara sum
                            addXTilde!(Workspace, is, it, iu, nw,T) # add to XTilde-type bubble functions
                        end
                    end
                end
            end
        end
    end
end

@inline function mixedFrequencies(ns, nt, nu, nwpr)
    nw1 = Int((ns + nt + nu - 1) / 2)
    nw2 = Int((ns - nt - nu - 1) / 2)
    nw3 = Int((-ns + nt - nu - 1) / 2)
    nw4 = Int((-ns - nt + nu - 1) / 2)
    wpw1 = nwpr + nw1 + 1
    wpw2 = nwpr + nw2 + 1
    wmw3 = nwpr - nw3
    wmw4 = nwpr - nw4
    # @assert (ns + wmw3 +wmw4)%2 != 0 "error in freq"
    return wpw1, wpw2, wmw3, wmw4
end

"""
adds part of X functions in Matsubara sum at nwpr containing the site summation for a set of s t and u frequencies. This is the most numerically demanding part!
"""
function addX!(
    Workspace::PMFRGWorkspace,
    Rij::Integer,
    is::Integer,
    it::Integer,
    iu::Integer,
    PropsBuffers,
    VertexBuffer,
)
    (; X, Par) = Workspace
    (; Va_ki_kj, Vb_ki_kj, Vc_ki_kj) = VertexBuffer
    (; N, np_vec,lenIntw) = Par.NumericalParams
    Nsum = Par.System.Nsum[Rij]
    ns = np_vec[is]
    nt = np_vec[it]
    nu = np_vec[iu]
    
    Xa_sum = zero(eltype(X.a)) #Perform summation on this temp variable before writing to State array as Base.setindex! proved to be a bottleneck!
    Xb_sum = zero(eltype(X.b))
    Xc_sum = zero(eltype(X.c))
    
    for (iw,nwpr) in enumerate(-lenIntw:lenIntw-1) # Matsubara sum
        wpw1, wpw2, wmw3, wmw4 = mixedFrequencies(ns, nt, nu, nwpr)
        PropsBuff = @view PropsBuffers[1:Nsum,iw]

        Va12 = getBufferView_ki(Va_ki_kj, ns, wpw1, wpw2, N, Nsum)
        Vb12 = getBufferView_ki(Vb_ki_kj, ns, wpw1, wpw2, N, Nsum)
        Vc12 = getBufferView_ki(Vc_ki_kj, ns, wpw1, wpw2, N, Nsum)

        Va34 = getBufferView_kj(Va_ki_kj, ns, wmw3, wmw4, N, Nsum)
        Vb34 = getBufferView_kj(Vb_ki_kj, ns, wmw3, wmw4, N, Nsum)
        Vc34 = getBufferView_kj(Vc_ki_kj, ns, wmw3, wmw4, N, Nsum)

        Vc21 = getBufferView_ki(Vc_ki_kj, ns, wpw2, wpw1, N, Nsum)
        Vc43 = getBufferView_kj(Vc_ki_kj, ns, wmw4, wmw3, N, Nsum)

        @turbo unroll = 1 for k in eachindex(Va12,Vb12,Vc12,Va34,Vb34,Vc34,Vc21,Vc43,PropsBuff)
            #loop over all Nsum summation elements defined in geometry. This inner loop is responsible for most of the computational effort! 
            Xa_sum += (+Va12[k] * Va34[k] + Vb12[k] * Vb34[k] * 2) * PropsBuff[k]

            Xb_sum +=
                (+Va12[k] * Vb34[k] + Vb12[k] * Va34[k] + Vb12[k] * Vb34[k]) * PropsBuff[k]

            Xc_sum += (+Vc12[k] * Vc34[k] + Vc21[k] * Vc43[k]) * PropsBuff[k]
        end
    end
    X.a[Rij, is, it, iu] = Xa_sum
    X.b[Rij, is, it, iu] = Xb_sum
    X.c[Rij, is, it, iu] = Xc_sum
end

function addXTilde!(
    Workspace::PMFRGWorkspace,
    is::Integer,
    it::Integer,
    iu::Integer,
    nwpr::Integer,
    T::AbstractFloat,
)

    (; State, X, Par) = Workspace
    (; N, np_vec) = Par.NumericalParams
    (; Npairs, invpairs, PairTypes, OnsitePairs) = Par.System

    @inline iG(x, nw) = iG_(Workspace.State.γ, x, T, nw)
    @inline iSKat(x, nw) = iSKat_(Workspace.State.γ, Workspace.Deriv.γ, x, T, nw)

    ns = np_vec[is]
    nt = np_vec[it]
    nu = np_vec[iu]
    
    @inline Props(i,j) = iSKat(i, nwpr) * iG(j, nwpr + ns)

    @inline Va_(Rij, s, t, u) = V_(State.Γ.a, Rij, s, t, u, invpairs[Rij], N)
    @inline Vb_(Rij, s, t, u) = V_(State.Γ.b, Rij, s, t, u, invpairs[Rij], N)
    @inline Vc_(Rij, s, t, u) = V_(State.Γ.c, Rij, s, t, u, invpairs[Rij], N)

    wpw1, wpw2, wmw3, wmw4 = mixedFrequencies(ns, nt, nu, nwpr)

    #Xtilde only defined for nonlocal pairs Rij != Rii
    for Rij = 1:Npairs
        Rij in OnsitePairs && continue
        #loop over all left hand side inequivalent pairs Rij
        Rji = invpairs[Rij] # store pair corresponding to Rji (easiest case: Rji = Rij)
        (; xi, xj) = PairTypes[Rij]

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

        X.Ta[Rij, is, it, iu] += (
            (+Va21 * Va43 + 2 * Vc21 * Vc43) * Props(xi, xj) +
            (Va12 * Va34 + 2 * Vc12 * Vc34) * Props(xj, xi)
        )

        X.Tb[Rij, is, it, iu] += (
            (+Va21 * Vc43 + Vc21 * Vc43 + Vc21 * Va43) * Props(xi, xj) +
            (Va12 * Vc34 + Vc12 * Vc34 + Vc12 * Va34) * Props(xj, xi)
        )
        Vb12 = Vb_(Rji, wpw1, wpw2, ns)
        Vb21 = Vb_(Rij, wpw2, wpw1, ns)
        Vb34 = Vb_(Rji, wmw3, wmw4, ns)
        Vb43 = Vb_(Rij, wmw4, wmw3, ns)

        Vc12 = Vc_(Rji, wpw1, wpw2, ns)
        Vc21 = Vc_(Rij, wpw2, wpw1, ns)
        Vc34 = Vc_(Rji, wmw3, wmw4, ns)
        Vc43 = Vc_(Rij, wmw4, wmw3, ns)


        X.Tc[Rij, is, it, iu] += (
            (+Vb21 * Vb43 + Vc21 * Vc43) * Props(xi, xj) +
            (Vb12 * Vb34 + Vc12 * Vc34) * Props(xj, xi)
        )
    end
end

"""Use symmetries and identities to compute the rest of bubble functions"""
function symmetrizeBubble!(X::BubbleType, Par::PMFRGParams)
    N = Par.NumericalParams.N
    (; Npairs, OnsitePairs) = Par.System
    usesymmetry = Par.Options.usesymmetry
    # use the u <--> t symmetry
    if (usesymmetry)
        Threads.@threads for it = 1:N
            for iu = it+1:N, is = 1:N, Rij = 1:Npairs
                X.a[Rij, is, it, iu] = -X.a[Rij, is, iu, it]
                X.b[Rij, is, it, iu] = -X.b[Rij, is, iu, it]
                X.c[Rij, is, it, iu] =
                    (+X.a[Rij, is, it, iu] - X.b[Rij, is, it, iu] + X.c[Rij, is, iu, it])
            end
        end
    end
    #local definitions of X.Tilde vertices
    Threads.@threads for iu = 1:N
        for it = 1:N, is = 1:N, R in OnsitePairs
            X.Ta[R, is, it, iu] = X.a[R, is, it, iu]
            X.Tb[R, is, it, iu] = X.b[R, is, it, iu]
            X.Tc[R, is, it, iu] = X.c[R, is, it, iu]
            X.Td[R, is, it, iu] = -X.c[R, is, iu, it]
        end
    end
    @tturbo X.Td .= X.Ta .- X.Tb .- X.Tc
end

function addToVertexFromBubble!(Γ::VertexType, X::BubbleType)
    Threads.@threads for iu in axes(Γ.a, 4)
        for it in axes(Γ.a, 3), is in axes(Γ.a, 2), Rij in axes(Γ.a, 1)
            Γ.a[Rij, is, it, iu] +=
                X.a[Rij, is, it, iu] - X.Ta[Rij, it, is, iu] + X.Ta[Rij, iu, is, it]
            Γ.b[Rij, is, it, iu] +=
                X.b[Rij, is, it, iu] - X.Tc[Rij, it, is, iu] + X.Tc[Rij, iu, is, it]
            Γ.c[Rij, is, it, iu] +=
                X.c[Rij, is, it, iu] - X.Tb[Rij, it, is, iu] + X.Td[Rij, iu, is, it]
        end
    end
    return Γ
end


function symmetrizeVertex!(Γ::VertexType, Par)
    N = Par.NumericalParams.N
    Threads.@threads for iu = 1:N
        for it = 1:N, is = 1:N, R in Par.System.OnsitePairs
            Γ.c[R, is, it, iu] = -Γ.b[R, it, is, iu]
        end
    end
end

##
getChi(State::ArrayPartition, T::Real, Par::PMFRGParams, Numax) =
    getChi(State.x[2], State.x[5], T, Par, Numax)
getChi(State::ArrayPartition, T::Real, Par::PMFRGParams) =
    getChi(State.x[2], State.x[5], T, Par)

function getChi(gamma::AbstractArray, Γc::AbstractArray, T::Real, Par::PMFRGParams, Numax)
    (; N, lenIntw_acc, np_vec) = Par.NumericalParams
    (; Npairs, invpairs, PairTypes, OnsitePairs) = Par.System

    @inline iG(x, w) = iG_(gamma, x, T, w)
    @inline Vc_(Rij, s, t, u) = V_(Γc, Rij, s, t, u, invpairs[Rij], N)

    Chi = zeros(_getFloatType(Par), Npairs, N)

    @inbounds Threads.@threads for Rij = 1:Npairs
        (; xi, xj) = PairTypes[Rij]
        for i_nu = 1:Numax
            n_nu = np_vec[i_nu]

            for nK = -lenIntw_acc:lenIntw_acc-1
                if Rij in OnsitePairs
                    Chi[Rij, i_nu] += iG(xi, nK) * iG(xi, nK + n_nu)
                end
                for nK2 = -lenIntw_acc:lenIntw_acc-1
                    npwpw2 = n_nu + nK + nK2 + 1
                    wmw2 = nK - nK2
                    #use that Vc_0 is calculated from Vb
                    GGGG = iG(xi, nK) * iG(xi, nK + n_nu) * iG(xj, nK2) * iG(xj, nK2 + n_nu)
                    Chi[Rij, i_nu] += GGGG * Vc_(Rij, n_nu, npwpw2, wmw2)
                end
            end
        end
    end
    return (Chi)
end

function getChi(gamma::AbstractArray, Γc::AbstractArray, T::Real, Par::PMFRGParams)
    (; N, lenIntw_acc) = Par.NumericalParams
    (; Npairs, invpairs, PairTypes, OnsitePairs) = Par.System

    @inline iG(x, w) = iG_(gamma, x, T, w)
    @inline Vc_(Rij, s, t, u) = V_(Γc, Rij, s, t, u, invpairs[Rij], N)

    Chi = zeros(_getFloatType(Par), Npairs)

    @inbounds Threads.@threads for Rij = 1:Npairs
        (; xi, xj) = PairTypes[Rij]
        for nK = -lenIntw_acc:lenIntw_acc-1
            if Rij in OnsitePairs
                Chi[Rij] += iG(xi, nK)^2
            end
            for nK2 = -lenIntw_acc:lenIntw_acc-1
                npwpw2 = nK + nK2 + 1
                wmw2 = nK - nK2
                #use that Vc_0 is calculated from Vb
                GGGG = iG(xi, nK)^2 * iG(xj, nK2)^2
                Chi[Rij] += GGGG * Vc_(Rij, 0, npwpw2, wmw2)
            end
        end
    end
    return (Chi)
end
