function test_BubbleSymmetries(B::BubbleType,BTranspose::BubbleType = B;kwargs...)
    test_tu_symmetry_ab(B.a,"Ba";kwargs...)
    test_tu_symmetry_ab(B.b,"Bb";kwargs...)
    test_tu_symmetry_c(B.a,B.b,B.c,"B";kwargs...)
end

# function Base.isequal(x::T,y::T) where {T <: Union{BubbleType,VertexType}}
#     for f in fieldnames(T)
#         for (i,j) in zip(getfield(x,f),getfield(y,f))
#             i == j || return false
#         end
#     end
#     return true
# end

function test_BubbleSymmetries(Workspace::ParquetWorkspace;kwargs...)
    @testset "test B0 Bubble" begin test_BubbleSymmetries(Workspace.B0) end
    @testset "test BX Bubble" begin test_BubbleSymmetries(Workspace.BX) end
end

# test_BubbleSymmetries(Par::ParquetParams = Params(getPolymer(2),Parquet(),N=24,lenIntw = 120);kwargs...) = test_BubbleSymmetries(test_computeBubbles(Par);kwargs...)

# Todo test that B0 bubble gives the same result as normal bubble with B0 inserted

function test_computeBubbles(Par)
    Workspace = SetupParquet(Par)
    @unpack State,Γ0,X,B0,BX,Par,Buffer = Workspace
    getProp! = constructPropagatorFunction(Workspace,1.)

    computeLeft2PartBubble!(B0,Γ0,Γ0,State.Γ,getProp!,Par,Buffer)
    computeLeft2PartBubble!(BX,X,X,State.Γ,getProp!,Par,Buffer)
    return Workspace
end

# function test_BareBubbles(Par::ParquetParams = Params(getPolymer(2),Parquet(),N=24,lenIntw = 120);tol = 1e-14)
#     WS = test_computeBubbles(Par)
#     @testset "B0 == BX(Γ^0,Γ)" begin @test dist(WS.B0,WS.BX) ≈ 0. atol = tol end
# end

function test_BareBubbles(;tol = 1e-14)
    
    Sol,Obs,Par = test_runDimerParquet()
    getProp! = constructPropagatorFunction(Sol,1.)
    
    @unpack State,Γ0,X,B0,BX,Par,Buffer = Sol

    Γinit = constructBubbleFromVertex(BareVertex_Freq(Par))
    computeLeft2PartBubble!(B0,Γinit,Γinit,State.Γ,getProp!,Par,Buffer)
    computeLeft2PartBubble!(BX,X,X,State.Γ,getProp!,Par,Buffer)
    @testset "B0 == BX(Γ^0,Γ)" begin @test dist(Sol.B0,Sol.BX) ≈ 0. atol = tol end
end