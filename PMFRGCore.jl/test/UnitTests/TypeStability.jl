
"""gets empty Workspace."""
function getEmptyWorkspace(Par::PMFRGCore.OneLoopParams)
    State = PMFRGCore.InitializeState(Par)
    Deriv = similar(State)
    setup = PMFRGCore.AllocateSetup(Par)
    (X, Buffs, Par) = setup
    PMFRGCore.OneLoopWorkspace(Deriv, State, X, Buffs, Par)
end

convertToSMatrix(M::PMFRGCore.MMatrix{N,N,T,NN}) where {N,T,NN} =
    PMFRGCore.SMatrix{N,N,T,NN}(M)

"""Checks whether vertex derivative has type instabilities"""
function test_OneLoopAllocations(Par = Params(getPolymer(2)))
    WS = getEmptyWorkspace(Par)
    # eval("@code_warntype PMFRGCore.getXBubble!(WS,0.)")
    # @code_warntype PMFRGCore.getXBubble!(WS,0.)
    SBuff = take!(WS.Buffer.Props)
    VBuff = take!(WS.Buffer.Vertex)
    SBuff1 = convertToSMatrix(SBuff)

    SBuffMatrix = Matrix(SBuff1)
    test_OneLoopAllocations(Par, WS, SBuff1, VBuff)
    test_OneLoopAllocations(Par, WS, SBuffMatrix, VBuff)

    put!(WS.Buffer.Props, SBuff)
    put!(WS.Buffer.Vertex, VBuff)
end
function test_OneLoopAllocations(Par, WS, SProps, VBuff)
    PMFRGCore.addX!(WS.X, WS.State, WS.Par, 1, 1, 2, 2, SProps, VBuff) # compile functions
    PMFRGCore.addXTilde!(WS.X, WS.State, WS.Par, 1, 1, 2, 2, SProps) # compile functions

    aX = @allocated PMFRGCore.addX!(WS.X, WS.State, WS.Par, 1, 1, 2, 2, SProps, VBuff)
    aXT = @allocated PMFRGCore.addXTilde!(WS.X, WS.State, WS.Par, 1, 1, 2, 2, SProps)
    test_TypeStability(Par, WS, aX, aXT)
end

function test_TypeStability(Par, WS, aX, aXT)
    @testset verbose = true "Allocations and Type stability" begin
        @testset "Type concreteness of Parameters " begin
            @test isconcretetype(typeof(Par))
        end
        @testset "Type concreteness of Workspace " begin
            @test isconcretetype(typeof(WS))
        end
        @testset "zero allocations of addX! " begin
            @test aX == 0
        end
        @testset "zero allocations of add XTilde! " begin
            @test aXT == 0
        end
    end
end

function getEmptyWorkspace(Par::PMFRGCore.TwoLoopParams)
    State = PMFRGCore.InitializeState(Par)
    Deriv = similar(State)
    setup = PMFRGCore.AllocateSetup(Par)
    (X, Y, Buffs, Par) = setup
    PMFRGCore.TwoLoopWorkspace(Deriv, State, X, Y, Buffs, Par)
end
function test_TwoLoopAllocations(Par = Params(getPolymer(2), TwoLoop()))
    WS = getEmptyWorkspace(Par)
    PropB_0 = take!(WS.Buffer.Props)

    VB = take!(WS.Buffer.Vertex)
    XB = take!(WS.Buffer.X)

    PropB = convertToSMatrix(PropB_0)
    PropBMatrix = Matrix(PropB)

    test_TwoLoopAllocations(Par, WS, PropB, VB, XB)
    test_TwoLoopAllocations(Par, WS, PropBMatrix, VB, XB)

    put!(WS.Buffer.Vertex, VB)
    put!(WS.Buffer.X, XB)
    put!(WS.Buffer.Props, PropB_0)
end
function test_TwoLoopAllocations(Par, WS, PropB, VB, XB)
    PMFRGCore.addBL!(WS.Y, WS.X, WS.X, WS.State.Γ, 1, 1, 2, 2, Par, PropB, VB, XB) #compile functions
    PMFRGCore.addBR!(WS.Y, WS.X, WS.X, WS.State.Γ, 1, 1, 2, 2, Par, PropB, VB, XB) #compile functions
    PMFRGCore.addBLTilde!(WS.Y, WS.X, WS.X, WS.State.Γ, 1, 1, 2, 2, Par, PropB) #compile functions
    PMFRGCore.addBRTilde!(WS.Y, WS.X, WS.X, WS.State.Γ, 1, 1, 2, 2, Par, PropB) #compile functions

    aX = @allocated PMFRGCore.addBL!(
        WS.Y,
        WS.X,
        WS.X,
        WS.State.Γ,
        1,
        1,
        2,
        2,
        Par,
        PropB,
        VB,
        XB,
    )
    aX += @allocated PMFRGCore.addBR!(
        WS.Y,
        WS.X,
        WS.X,
        WS.State.Γ,
        1,
        1,
        2,
        2,
        Par,
        PropB,
        VB,
        XB,
    )
    aXT = @allocated PMFRGCore.addBLTilde!(
        WS.Y,
        WS.X,
        WS.X,
        WS.State.Γ,
        1,
        1,
        2,
        2,
        Par,
        PropB,
    )
    aXT += @allocated PMFRGCore.addBRTilde!(
        WS.Y,
        WS.X,
        WS.X,
        WS.State.Γ,
        1,
        1,
        2,
        2,
        Par,
        PropB,
    )

    test_TypeStability(Par, WS, aX, aXT)

end
# Todo: try to write unit tests using @code_warntype  

function capture_stdout(f)
    stdout_orig = stdout
    (rd, wr) = redirect_stdout()
    f()
    close(wr)
    redirect_stdout(stdout_orig)
    read(rd, String)
end
