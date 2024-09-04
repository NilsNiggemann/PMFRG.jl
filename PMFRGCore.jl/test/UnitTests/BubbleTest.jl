function test_BubbleSymmetries(
    B::PMFRGCore.BubbleType,
    BTranspose::PMFRGCore.BubbleType = B;
    kwargs...,
)
    test_tu_symmetry_ab(B.a, "Ba"; kwargs...)
    test_tu_symmetry_ab(B.b, "Bb"; kwargs...)
    test_tu_symmetry_c(B.a, B.b, B.c, "B"; kwargs...)
end
function test_BubbleSymmetries(Workspace::PMFRGCore.ParquetWorkspace; kwargs...)
    @testset "test B0 Bubble" begin
        test_BubbleSymmetries(Workspace.B0)
    end
    @testset "test BX Bubble" begin
        test_BubbleSymmetries(Workspace.BX)
    end
end

# test_BubbleSymmetries(Par::ParquetParams = Params(getPolymer(2),Parquet(),N=24,lenIntw = 120);kwargs...) = test_BubbleSymmetries(test_computeBubbles(Par);kwargs...)

# Todo test that B0 bubble gives the same result as normal bubble with B0 inserted

function test_computeBubbles(Par)
    Workspace = PMFRGCore.SetupParquet(Par)
    (; State, Γ0, X, B0, BX, Par, Buffer = Workspace)
    PMFRGCore.getProp! = constructPropagatorFunction(Workspace, 1.0)

    PMFRGCore.addToLeft2PartBubble!(B0, Γ0, Γ0, State.Γ, getProp!, Par, Buffer)
    PMFRGCore.addToLeft2PartBubble!(BX, X, X, State.Γ, getProp!, Par, Buffer)
    return Workspace
end

# function test_BareBubbles(Par::ParquetParams = Params(getPolymer(2),Parquet(),N=24,lenIntw = 120);tol = 1e-14)
#     WS = test_computeBubbles(Par)
#     @testset "B0 == BX(Γ^0,Γ)" begin @test dist(WS.B0,WS.BX) ≈ 0. atol = tol end
# end

function test_BareBubbles(; kwargs...)

    Sol, Obs, Par = test_runDimerParquet()
    test_BareBubbles(Sol; kwargs...)
end

function test_BareBubbles(Workspace::PMFRGCore.ParquetWorkspace; tol = 1e-14)

    getProp! = PMFRGCore.constructPropagatorFunction(Workspace, 1.0)

    (; State, Γ0, X, B0, BX, Par, Buffer) = Workspace

    B0_new = deepcopy(B0) # Allocate new memory to avoid mutating Workspace
    BX_new = deepcopy(BX)

    Γinit = PMFRGCore.constructBubbleFromVertex(PMFRGCore.BareVertex_Freq(Par))

    PMFRGCore.computeLeft2PartBubble!(B0_new, Γ0, Γ0, State.Γ, getProp!, Par, Buffer)
    PMFRGCore.computeLeft2PartBubble!(BX_new, Γinit, Γinit, State.Γ, getProp!, Par, Buffer)

    @testset "B0 == BX(Γ^0,Γ)" begin
        @test PMFRGCore.dist(B0_new, BX_new) ≈ 0.0 atol = tol
    end
end
