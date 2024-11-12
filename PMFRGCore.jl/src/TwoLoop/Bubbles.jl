
"""
adds to ResultBubble given the vertex as well as a bubble inserted on the left. Assumes that vertices and bubbles are given pre-computed in VertexBuffer. Providing the appropriate functions BufferFiller! and SiteIndex allows for computation of right insertions instead, since these are shown to only swap site and frequency indices.
"""

function addBL!(
    B::BubbleType{T},
    XL::BubbleType,
    XR::BubbleType,
    Γ::VertexType,
    is::Integer,
    it::Integer,
    iu::Integer,
    nwpr::Integer,
    Par::PMFRGParams,
    Props,
    VBuffer::VertexBufferType,
    XBuffer::BubbleBufferType,
    BufferFiller!::Func1 = fillBufferL!,
    SiteIndex::Func2 = identity,
) where {T,Func1<:Function,Func2<:Function}

    (; N, np_vec) = Par.NumericalParams
    (; Npairs, Nsum, siteSum, invpairs) = Par.System

    BufferFiller!(VBuffer, XBuffer, XL, XR, Γ, is, it, iu, nwpr, Par)

    (; Va34, Vb34, Vc34, Vc43) = VBuffer
    (; XTa21, XTb21, XTc21, XTd21) = XBuffer

    S_ki = siteSum.ki
    S_kj = siteSum.kj
    S_xk = siteSum.xk
    S_m = siteSum.m

    @inbounds for Rij = 1:Npairs
        #loop over all left hand side inequivalent pairs Rij
        Ba_sum = zero(T) #Perform summation on this temp variable before writing to State array as Base.setindex! proved to be a bottleneck!
        Bb_sum = zero(T)
        Bc_sum = zero(T)
        Site = SiteIndex(Rij)
        @turbo unroll = 1 for k_spl = 1:Nsum[Rij]
            #loop over all Nsum summation elements defined in geometry. This inner loop is responsible for most of the computational effort! 
            ki, kj, m, xk =
                S_ki[k_spl, Site], S_kj[k_spl, Site], S_m[k_spl, Site], S_xk[k_spl, Site]

            Ptm = Props[xk, xk] * m

            Ba_sum += (Va34[kj] * XTa21[ki] + 2 * Vb34[kj] * XTc21[ki]) * Ptm

            Bb_sum +=
                (Vb34[kj] * XTa21[ki] + Va34[kj] * XTc21[ki] + Vb34[kj] * XTc21[ki]) * Ptm


            Bc_sum += (-Vc43[kj] * XTb21[ki] + Vc34[kj] * XTd21[ki]) * Ptm
        end
        B.a[Rij, is, it, iu] += Ba_sum
        B.b[Rij, is, it, iu] += Bb_sum
        B.c[Rij, is, it, iu] += Bc_sum
    end
    return
end

function addBLTilde!(
    B::BubbleType,
    XL::BubbleType,
    XR::BubbleType,
    Γ::VertexType,
    is::Integer,
    it::Integer,
    iu::Integer,
    nwpr::Integer,
    Par::PMFRGParams,
    Props,
    FreqSwap::Func = identity,
) where {Func<:Function}

    (; Npairs, invpairs, PairTypes, OnsitePairs) = Par.System
    (; N, np_vec) = Par.NumericalParams

    @inline Va_(Rij, s, t, u) = V_(Γ.a, Rij, s, t, u, invpairs[Rij], N)
    @inline Vb_(Rij, s, t, u) = V_(Γ.b, Rij, s, t, u, invpairs[Rij], N)
    @inline Vc_(Rij, s, t, u) = V_(Γ.c, Rij, s, t, u, invpairs[Rij], N)

    @inline XLa_(Rij, s, t, u) = X_(XL.a, XR.a, Rij, s, t, u, invpairs[Rij], N)
    @inline XLb_(Rij, s, t, u) = X_(XL.b, XR.b, Rij, s, t, u, invpairs[Rij], N)
    @inline XLc_(Rij, s, t, u) = X_(XL.c, XR.c, Rij, s, t, u, invpairs[Rij], N)

    @inline XRTa_(Rij, s, t, u) = XT_(XR.Ta, XL.Ta, Rij, s, t, u, invpairs[Rij], N)
    @inline XRTb_(Rij, s, t, u) = XT_(XR.Tb, XL.Tb, Rij, s, t, u, invpairs[Rij], N)
    @inline XRTc_(Rij, s, t, u) = XT_(XR.Tc, XL.Tc, Rij, s, t, u, invpairs[Rij], N)
    @inline XRTd_(Rij, s, t, u) = XT_(XR.Td, XL.Td, Rij, s, t, u, invpairs[Rij], N)

    ns = np_vec[is]
    nt = np_vec[it]
    nu = np_vec[iu]
    wpw1, wpw2, wmw3, wmw4 = FreqSwap(mixedFrequencies(ns, nt, nu, nwpr))  #swap frequencies for right bubble
    #Btilde only defined for nonlocal pairs Rij != Rii
    for Rij = 1:Npairs
        Rij in OnsitePairs && continue
        #loop over all left hand side inequivalent pairs Rij
        Rji = invpairs[Rij] # store pair corresponding to Rji (easiest case: Rji = Rij)

        (; xi, xj) = PairTypes[Rij]

        B.Ta[Rij, is, it, iu] +=
            Props[xj, xi] * (
                Va_(Rij, wmw4, ns, wmw3) * XLa_(Rij, wpw2, ns, wpw1) +
                2 * Vc_(Rij, wmw4, ns, wmw3) * XLc_(Rij, wpw2, ns, wpw1)
            ) +
            Props[xi, xj] * (
                Va_(Rji, wmw3, ns, wmw4) * XRTa_(Rji, wpw2, wpw1, ns) +
                2 * Vc_(Rji, wmw3, ns, wmw4) * XRTd_(Rji, wpw2, wpw1, ns)
            )

        B.Tb[Rij, is, it, iu] +=
            Props[xj, xi] * (
                Va_(Rij, wmw4, ns, wmw3) * XLc_(Rij, wpw2, ns, wpw1) +
                Vc_(Rij, wmw4, ns, wmw3) *
                (XLa_(Rij, wpw2, ns, wpw1) + XLc_(Rij, wpw2, ns, wpw1))
            ) +
            Props[xi, xj] * (
                Va_(Rji, wmw3, ns, wmw4) * XRTd_(Rji, wpw2, wpw1, ns) +
                Vc_(Rji, wmw3, ns, wmw4) *
                (XRTa_(Rji, wpw2, wpw1, ns) + XRTd_(Rji, wpw2, wpw1, ns))
            )

        B.Tc[Rij, is, it, iu] +=
            Props[xj, xi] * (
                Vb_(Rij, wmw4, ns, wmw3) * XLb_(Rij, wpw2, ns, wpw1) +
                Vc_(Rij, wmw4, wmw3, ns) * XLc_(Rij, wpw2, wpw1, ns)
            ) +
            Props[xi, xj] * (
                -Vc_(Rji, wmw3, wmw4, ns) * XRTb_(Rji, wpw2, wpw1, ns) +
                Vb_(Rji, wmw3, ns, wmw4) * XRTc_(Rji, wpw2, wpw1, ns)
            )
    end
end

function addBL!(
    B::BubbleType{T},
    XL::BubbleType,
    XR::BubbleType,
    Γ::VertexType,
    is::Integer,
    it::Integer,
    iu::Integer,
    nwpr::Integer,
    Par::PMFRGParams,
    Props::SingleElementMatrix,
    VBuffer::VertexBufferType,
    XBuffer::BubbleBufferType,
    BufferFiller!::Func1 = fillBufferL!,
    SiteIndex::Func2 = identity,
) where {T,Func1<:Function,Func2<:Function}

    (; Npairs, Nsum, siteSum) = Par.System

    BufferFiller!(VBuffer, XBuffer, XL, XR, Γ, is, it, iu, nwpr, Par)

    (; Va34, Vb34, Vc34, Vc43) = VBuffer
    (; XTa21, XTb21, XTc21, XTd21) = XBuffer

    S_ki = siteSum.ki
    S_kj = siteSum.kj
    S_m = siteSum.m
    Prop = only(Props)

    @inbounds for Rij = 1:Npairs
        #loop over all left hand side inequivalent pairs Rij
        Ba_sum = zero(T) #Perform summation on this temp variable before writing to State array as Base.setindex! proved to be a bottleneck!
        Bb_sum = zero(T)
        Bc_sum = zero(T)
        Site = SiteIndex(Rij)
        @turbo unroll = 1 for k_spl = 1:Nsum[Rij]
            #loop over all Nsum summation elements defined in geometry. This inner loop is responsible for most of the computational effort! 
            ki, kj, m = S_ki[k_spl, Site], S_kj[k_spl, Site], S_m[k_spl, Site]

            mConv = convert(T, m)
            Ba_sum += (Va34[kj] * XTa21[ki] + 2 * Vb34[kj] * XTc21[ki]) * mConv

            Bb_sum +=
                (Vb34[kj] * XTa21[ki] + Va34[kj] * XTc21[ki] + Vb34[kj] * XTc21[ki]) * mConv


            Bc_sum += (-Vc43[kj] * XTb21[ki] + Vc34[kj] * XTd21[ki]) * mConv
        end
        B.a[Rij, is, it, iu] += Ba_sum * Prop
        B.b[Rij, is, it, iu] += Bb_sum * Prop
        B.c[Rij, is, it, iu] += Bc_sum * Prop
    end
    return
end

function addBLTilde!(
    B::BubbleType,
    XL::BubbleType,
    XR::BubbleType,
    Γ::VertexType,
    is::Integer,
    it::Integer,
    iu::Integer,
    nwpr::Integer,
    Par::PMFRGParams,
    Props::SingleElementMatrix,
    FreqSwap::Func = identity,
) where {Func<:Function}

    (; Npairs, invpairs, PairTypes, OnsitePairs) = Par.System
    (; N, np_vec) = Par.NumericalParams

    @inline Va_(Rij, s, t, u) = V_(Γ.a, Rij, s, t, u, invpairs[Rij], N)
    @inline Vb_(Rij, s, t, u) = V_(Γ.b, Rij, s, t, u, invpairs[Rij], N)
    @inline Vc_(Rij, s, t, u) = V_(Γ.c, Rij, s, t, u, invpairs[Rij], N)

    @inline XLa_(Rij, s, t, u) = X_(XL.a, XR.a, Rij, s, t, u, invpairs[Rij], N)
    @inline XLb_(Rij, s, t, u) = X_(XL.b, XR.b, Rij, s, t, u, invpairs[Rij], N)
    @inline XLc_(Rij, s, t, u) = X_(XL.c, XR.c, Rij, s, t, u, invpairs[Rij], N)

    @inline XRTa_(Rij, s, t, u) = XT_(XR.Ta, XL.Ta, Rij, s, t, u, invpairs[Rij], N)
    @inline XRTb_(Rij, s, t, u) = XT_(XR.Tb, XL.Tb, Rij, s, t, u, invpairs[Rij], N)
    @inline XRTc_(Rij, s, t, u) = XT_(XR.Tc, XL.Tc, Rij, s, t, u, invpairs[Rij], N)
    @inline XRTd_(Rij, s, t, u) = XT_(XR.Td, XL.Td, Rij, s, t, u, invpairs[Rij], N)

    ns = np_vec[is]
    nt = np_vec[it]
    nu = np_vec[iu]
    wpw1, wpw2, wmw3, wmw4 = FreqSwap(mixedFrequencies(ns, nt, nu, nwpr))  #swap frequencies for right bubble
    Prop = only(Props)
    #Btilde only defined for nonlocal pairs Rij != Rii
    for Rij = 1:Npairs
        Rij in OnsitePairs && continue
        #loop over all left hand side inequivalent pairs Rij
        Rji = invpairs[Rij] # store pair corresponding to Rji (easiest case: Rji = Rij)

        B.Ta[Rij, is, it, iu] +=
            Prop * (
                Va_(Rij, wmw4, ns, wmw3) * XLa_(Rij, wpw2, ns, wpw1) +
                2 * Vc_(Rij, wmw4, ns, wmw3) * XLc_(Rij, wpw2, ns, wpw1) +
                Va_(Rji, wmw3, ns, wmw4) * XRTa_(Rji, wpw2, wpw1, ns) +
                2 * Vc_(Rji, wmw3, ns, wmw4) * XRTd_(Rji, wpw2, wpw1, ns)
            )

        B.Tb[Rij, is, it, iu] +=
            Prop * (
                Va_(Rij, wmw4, ns, wmw3) * XLc_(Rij, wpw2, ns, wpw1) +
                Vc_(Rij, wmw4, ns, wmw3) *
                (XLa_(Rij, wpw2, ns, wpw1) + XLc_(Rij, wpw2, ns, wpw1)) +
                Va_(Rji, wmw3, ns, wmw4) * XRTd_(Rji, wpw2, wpw1, ns) +
                Vc_(Rji, wmw3, ns, wmw4) *
                (XRTa_(Rji, wpw2, wpw1, ns) + XRTd_(Rji, wpw2, wpw1, ns))
            )

        B.Tc[Rij, is, it, iu] +=
            Prop * (
                Vb_(Rij, wmw4, ns, wmw3) * XLb_(Rij, wpw2, ns, wpw1) +
                Vc_(Rij, wmw4, wmw3, ns) * XLc_(Rij, wpw2, wpw1, ns) -
                Vc_(Rji, wmw3, wmw4, ns) * XRTb_(Rji, wpw2, wpw1, ns) +
                Vb_(Rji, wmw3, ns, wmw4) * XRTc_(Rji, wpw2, wpw1, ns)
            )
    end
end

##
"""
adds to ResultBubble given the vertex as well as a bubble inserted on the left. Assumes that vertices and bubbles are given pre-computed in VertexBuffer.
"""
@inline function addBR!(
    B::BubbleType,
    XL::BubbleType,
    XR::BubbleType,
    Γ::VertexType,
    is::Integer,
    it::Integer,
    iu::Integer,
    nwpr::Integer,
    Par::PMFRGParams,
    Props,
    VBuffer::VertexBufferType,
    XBuffer::BubbleBufferType,
)
    addBL!(
        B::BubbleType,
        XL::BubbleType,
        XR::BubbleType,
        Γ::VertexType,
        is::Integer,
        it::Integer,
        iu::Integer,
        nwpr::Integer,
        Par::PMFRGParams,
        Props,
        VBuffer::VertexBufferType,
        XBuffer::BubbleBufferType,
        fillBufferR!,
        x -> Par.System.invpairs[x],
    )
end

@inline function addBRTilde!(
    B::BubbleType,
    XL::BubbleType,
    XR::BubbleType,
    Γ::VertexType,
    is::Integer,
    it::Integer,
    iu::Integer,
    nwpr::Integer,
    Par::PMFRGParams,
    Props,
)
    addBLTilde!(
        B::BubbleType,
        XL::BubbleType,
        XR::BubbleType,
        Γ::VertexType,
        is::Integer,
        it::Integer,
        iu::Integer,
        nwpr::Integer,
        Par::PMFRGParams,
        Props,
        freqs -> (freqs[3], freqs[4], freqs[1], freqs[2]),
    )
end
##

function addBL!(
    B::BubbleType{T},
    Γ0::BareVertexType,
    Γ::VertexType,
    is::Integer,
    it::Integer,
    iu::Integer,
    nwpr::Integer,
    Par::PMFRGParams,
    Props,
    Buffer::VertexBufferType,
) where {T}
    (; Npairs, Nsum, siteSum) = Par.System

    fillBufferL!(Buffer, Γ0, Γ, is, it, iu, nwpr, Par)

    (; Vc34, Vc43) = Buffer
    S_ki = siteSum.ki
    S_kj = siteSum.kj
    S_xk = siteSum.xk
    S_m = siteSum.m

    @inbounds for Rij = 1:Npairs
        Bc_sum = zero(T)
        @turbo unroll = 1 for k_spl = 1:Nsum[Rij]
            #loop over all Nsum summation elements defined in geometry. This inner loop is responsible for most of the computational effort! 
            ki, kj, m, xk =
                S_ki[k_spl, Rij], S_kj[k_spl, Rij], S_m[k_spl, Rij], S_xk[k_spl, Rij]

            Ptm = Props[xk, xk] * m

            Bc_sum += (Vc43[kj] + Vc34[kj]) * Ptm * Γ0.c[ki]
        end
        B.c[Rij, is, it, iu] += Bc_sum
    end
    return
end

@inline addBL!(
    B::BubbleType,
    Γ0L::BareVertexType,
    Γ0R::BareVertexType,
    Γ::VertexType,
    is::Integer,
    it::Integer,
    iu::Integer,
    nwpr::Integer,
    Par::PMFRGParams,
    Props,
    Buffer::VertexBufferType,
    XB::BubbleBufferType,
) = addBL!(B, Γ0L, Γ, is, it, iu, nwpr, Par, Props, Buffer)


function addBLTilde!(
    B::BubbleType,
    Γ0::BareVertexType,
    Γ::VertexType,
    is::Integer,
    it::Integer,
    iu::Integer,
    nwpr::Integer,
    Par::PMFRGParams,
    Props,
)

    (; Npairs, invpairs, PairTypes, OnsitePairs) = Par.System
    (; N, np_vec) = Par.NumericalParams


    @inline Va_(Rij, s, t, u) = V_(Γ.a, Rij, s, t, u, invpairs[Rij], N)
    @inline Vb_(Rij, s, t, u) = V_(Γ.b, Rij, s, t, u, invpairs[Rij], N)
    @inline Vc_(Rij, s, t, u) = V_(Γ.c, Rij, s, t, u, invpairs[Rij], N)

    ns = np_vec[is]
    nt = np_vec[it]
    nu = np_vec[iu]
    wpw1, wpw2, wmw3, wmw4 = mixedFrequencies(ns, nt, nu, nwpr)

    #Btilde only defined for nonlocal pairs Rij != Rii
    for Rij = 1:Npairs
        Rij in OnsitePairs && continue
        #loop over all left hand side inequivalent pairs Rij
        Rji = invpairs[Rij] # store pair corresponding to Rji (easiest case: Rji = Rij)
        (; xi, xj) = PairTypes[Rij]

        B.Ta[Rij, is, it, iu] +=
            Props[xj, xi] * 2 * Vc_(Rij, wmw4, ns, wmw3) * Γ0.c[Rij] +
            Props[xi, xj] * 2 * Vc_(Rji, wmw3, ns, wmw4) * Γ0.c[Rji]

        B.Tb[Rij, is, it, iu] +=
            Props[xj, xi] *
            (Va_(Rij, wmw4, ns, wmw3) + Vc_(Rij, wmw4, ns, wmw3)) *
            Γ0.c[Rij] +
            Props[xi, xj] *
            (Va_(Rji, wmw3, ns, wmw4) + Vc_(Rji, wmw3, ns, wmw4)) *
            Γ0.c[Rji]

        B.Tc[Rij, is, it, iu] +=
            Props[xj, xi] * Vc_(Rij, wmw4, wmw3, ns) * Γ0.c[Rij] +
            Props[xi, xj] * Vc_(Rji, wmw3, wmw4, ns) * Γ0.c[Rji]
    end
end
@inline addBLTilde!(
    B::BubbleType,
    Γ0L::BareVertexType,
    Γ0R::BareVertexType,
    Γ::VertexType,
    is::Integer,
    it::Integer,
    iu::Integer,
    nwpr::Integer,
    Par::PMFRGParams,
    Props,
) = addBLTilde!(B, Γ0L, Γ, is, it, iu, nwpr, Par, Props)

function compute1PartBubble_BS!(Dgamma::AbstractArray, Γ, Gamma0, Prop, Par)
    setZero!(Dgamma)

    (; T, Ngamma, lenIntw_acc, np_vec_gamma) = Par.NumericalParams
    (; siteSum, invpairs, Nsum, OnsitePairs) = Par.System

    @inline Γc_(Rij, s, t, u) = V_(Γ.c, Rij, s, t, u, invpairs[Rij], Par.NumericalParams.N) # cTilde corresponds to b type vertex!


    Threads.@threads for iw1 = 1:Ngamma
        nw1 = np_vec_gamma[iw1]
        for (x, Rx) in enumerate(OnsitePairs)
            jsum = 0.0
            for k_spl = 1:Nsum[Rx]
                (; m, ki, xk) = siteSum[k_spl, Rx]
                ik = invpairs[ki]
                fac = ifelse(ik in OnsitePairs, 1, 3)
                for nw = -lenIntw_acc:lenIntw_acc-1
                    for nwpr = -lenIntw_acc:lenIntw_acc-1

                        wpwpr = nw + nwpr + 1
                        w1mw = nw1 - nw
                        w1mwpr = nw1 - nwpr

                        jsum +=
                            fac *
                            4 *
                            Γc_(ik, wpwpr, w1mw, w1mwpr) *
                            Gamma0.c[ik] *
                            Prop(xk, nw) *
                            Prop(xk, nwpr) *
                            Prop(x, w1mw - nwpr - 1) *
                            m
                    end
                end
            end
            Dgamma[x, iw1] += -T^2 * jsum * 6^2 #6^2 since each Prop contains factor of 1/6
        end
    end
    return Dgamma
end
