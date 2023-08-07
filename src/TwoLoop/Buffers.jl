

@inline function X_(
    X::AbstractArray,
    XTransp::AbstractArray,
    Rj::Integer,
    ns::Integer,
    nt::Integer,
    nu::Integer,
    Rji::Integer,
    N::Integer,
)
    # @assert (ns+nt+nu) %2 != 0 "$ns + $nt +  $nu = $(ns+nt+nu)"
    ns, nt, nu, swapsites = convertFreqArgs(ns, nt, nu, N)
    Vertex = ifelse(swapsites, XTransp, X)
    Rj = ifelse(swapsites, Rji, Rj)
    return @inbounds Vertex[Rj, ns+1, nt+1, nu+1]
end

@inline function XT_(
    XT::AbstractArray,
    XTTransp::AbstractArray,
    Rj::Integer,
    ns::Integer,
    nt::Integer,
    nu::Integer,
    Rji::Integer,
    N::Integer,
)
    # @assert (ns+nt+nu) %2 != 0 "$ns + $nt +  $nu = $(ns+nt+nu)"
    Vertex = ifelse(nt * nu < 0, XTTransp, XT)
    ns, nt, nu, swapsites = convertFreqArgsXT(ns, nt, nu, N)
    Rj = ifelse(swapsites, Rji, Rj)
    return @inbounds Vertex[Rj, ns+1, nt+1, nu+1]
end

@inline function bufferXT_!(
    Cache,
    X::AbstractArray,
    XTransp::AbstractArray,
    ns::Integer,
    nt::Integer,
    nu::Integer,
    invpairs::AbstractArray,
    N,
)

    Vertex = ifelse(nt * nu < 0, XTransp, X)
    ns, nt, nu, swapsites = convertFreqArgsXT(ns, nt, nu, N)
    # @assert (ns+nt+nu) %2 != 0 "$ns + $nt +  $nu = $(ns+nt+nu)"

    is, it, iu = ns + 1, nt + 1, nu + 1
    @inbounds begin
        if swapsites
            @turbo unroll = 1 inline = true for R in eachindex(Cache, invpairs)
                Cache[R] = Vertex[invpairs[R], is, it, iu]
            end
        else
            @turbo unroll = 1 inline = true for R in eachindex(Cache, invpairs)
                Cache[R] = Vertex[R, is, it, iu]
            end
        end
    end
end

@inline function fillBufferL!(
    VBuffer::VertexBufferType,
    XBuffer::BubbleBufferType,
    XL::BubbleType,
    XR::BubbleType,
    Γ::VertexType,
    is::Integer,
    it::Integer,
    iu::Integer,
    nwpr::Integer,
    Par::PMFRGParams,
)
    np_vec = Par.NumericalParams.np_vec
    ns = np_vec[is]
    nt = np_vec[it]
    nu = np_vec[iu]
    wpw1, wpw2, wmw3, wmw4 = mixedFrequencies(ns, nt, nu, nwpr)

    fillBufferL_mixedFreq!(VBuffer, XBuffer, XL, XR, Γ, ns, wpw1, wpw2, wmw3, wmw4, Par)
end

@inline function fillBufferL_mixedFreq!(
    VBuffer::VertexBufferType,
    XBuffer::BubbleBufferType,
    XL::BubbleType,
    XR::BubbleType,
    Γ::VertexType,
    ns::Integer,
    wpw1::Integer,
    wpw2::Integer,
    wmw3::Integer,
    wmw4::Integer,
    Par::PMFRGParams,
)
    invpairs = Par.System.invpairs
    N = Par.NumericalParams.N


    (; Va34, Vb34, Vc34, Vc43) = VBuffer

    (; XTa21, XTb21, XTc21, XTd21) = XBuffer

    bufferV_!(Va34, Γ.a, ns, wmw3, wmw4, invpairs, N)
    bufferV_!(Vb34, Γ.b, ns, wmw3, wmw4, invpairs, N)
    bufferV_!(Vc34, Γ.c, ns, wmw3, wmw4, invpairs, N)

    bufferV_!(Vc43, Γ.c, ns, wmw4, wmw3, invpairs, N)

    bufferXT_!(XTa21, XR.Ta, XL.Ta, wpw2, ns, wpw1, invpairs, N)
    bufferXT_!(XTb21, XR.Tb, XL.Tb, wpw2, ns, wpw1, invpairs, N)
    bufferXT_!(XTc21, XR.Tc, XL.Tc, wpw2, ns, wpw1, invpairs, N)
    bufferXT_!(XTd21, XR.Td, XL.Td, wpw2, ns, wpw1, invpairs, N)
end


@inline function fillBufferR!(
    VBuffer::VertexBufferType,
    XBuffer::BubbleBufferType,
    XL::BubbleType,
    XR::BubbleType,
    Γ::VertexType,
    is::Integer,
    it::Integer,
    iu::Integer,
    nwpr::Integer,
    Par::PMFRGParams,
)

    np_vec = Par.NumericalParams.np_vec

    ns = np_vec[is]
    nt = np_vec[it]
    nu = np_vec[iu]
    wpw1, wpw2, wmw3, wmw4 = mixedFrequencies(ns, nt, nu, nwpr)
    wpw1, wpw2, wmw3, wmw4 = wmw3, wmw4, wpw1, wpw2 #swap frequencies for right bubble

    fillBufferL_mixedFreq!(VBuffer, XBuffer, XL, XR, Γ, ns, wpw1, wpw2, wmw3, wmw4, Par)
end

@inline function fillBufferL!(
    VBuffer::VertexBufferType,
    Γ0::BareVertexType,
    Γ::VertexType,
    is::Integer,
    it::Integer,
    iu::Integer,
    nwpr::Integer,
    Par::PMFRGParams,
)

    invpairs = Par.System.invpairs
    (; N, np_vec) = Par.NumericalParams

    (; Va34, Vb34, Vc34, Vc43) = VBuffer
    ns = np_vec[is]
    nt = np_vec[it]
    nu = np_vec[iu]
    wpw1, wpw2, wmw3, wmw4 = mixedFrequencies(ns, nt, nu, nwpr)

    bufferV_!(Va34, Γ.a, ns, wmw3, wmw4, invpairs, N)
    bufferV_!(Vb34, Γ.b, ns, wmw3, wmw4, invpairs, N)
    bufferV_!(Vc34, Γ.c, ns, wmw3, wmw4, invpairs, N)

    bufferV_!(Vc43, Γ.c, ns, wmw4, wmw3, invpairs, N)
end
