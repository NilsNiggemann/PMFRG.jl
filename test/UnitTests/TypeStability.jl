
"""gets empty Workspace."""
function getEmptyWorkspace(Par::PMFRG.OneLoopParams)
    State = PMFRG.InitializeState(Par)
    Deriv = similar(State)
    setup = PMFRG.AllocateSetup(Par)
    (X, Buffs, Par) = setup
    PMFRG.OneLoopWorkspace(Deriv, State, X, Buffs, Par)
end

convertToSMatrix(M::PMFRG.MMatrix{N,N,T,NN}) where {N,T,NN} = PMFRG.SMatrix{N,N,T,NN}(M)

"""Checks whether vertex derivative has type instabilities"""
function test_OneLoopAllocations(Par = Params(getPolymer(2)))
    WS = getEmptyWorkspace(Par)
    # eval("@code_warntype PMFRG.getXBubble!(WS,0.)")
    # @code_warntype PMFRG.getXBubble!(WS,0.)
    SBuff = take!(WS.Buffer.Props)
    VBuff = take!(WS.Buffer.Vertex)
    SBuff1 = convertToSMatrix(SBuff)
    test_OneLoopAllocations(Par, WS, SBuff1, VBuff)
    put!(WS.Buffer.Props, SBuff)
    put!(WS.Buffer.Vertex, VBuff)
end
function test_OneLoopAllocations(Par, WS, SProps, VBuff)
    PMFRG.addX!(WS.X, WS.State, WS.Par, 1, 1, 2, 2, SProps, VBuff) # compile functions
    PMFRG.addXTilde!(WS.X, WS.State, WS.Par, 1, 1, 2, 2, SProps) # compile functions

    aX = @allocated PMFRG.addX!(WS.X, WS.State, WS.Par, 1, 1, 2, 2, SProps, VBuff)
    aXT = @allocated PMFRG.addXTilde!(WS.X, WS.State, WS.Par, 1, 1, 2, 2, SProps)
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

function getEmptyWorkspace(Par::PMFRG.TwoLoopParams)
    State = PMFRG.InitializeState(Par)
    Deriv = similar(State)
    setup = PMFRG.AllocateSetup(Par)
    (X, Y, Buffs, Par) = setup
    PMFRG.TwoLoopWorkspace(Deriv, State, X, Y, Buffs, Par)
end
function test_TwoLoopAllocations(Par = Params(getPolymer(2), TwoLoop()))
    WS = getEmptyWorkspace(Par)
    PropB_0 = take!(WS.Buffer.Props)

    VB = take!(WS.Buffer.Vertex)
    XB = take!(WS.Buffer.X)

    PropB = convertToSMatrix(PropB_0)

    test_TwoLoopAllocations(Par, WS, PropB, VB, XB)
    put!(WS.Buffer.Vertex, VB)
    put!(WS.Buffer.X, XB)
    put!(WS.Buffer.Props, PropB_0)
end
function test_TwoLoopAllocations(Par, WS, PropB, VB, XB)
    PMFRG.addBL!(WS.Y, WS.X, WS.X, WS.State.Γ, 1, 1, 2, 2, Par, PropB, VB, XB) #compile functions
    PMFRG.addBR!(WS.Y, WS.X, WS.X, WS.State.Γ, 1, 1, 2, 2, Par, PropB, VB, XB) #compile functions
    PMFRG.addBLTilde!(WS.Y, WS.X, WS.X, WS.State.Γ, 1, 1, 2, 2, Par, PropB) #compile functions
    PMFRG.addBRTilde!(WS.Y, WS.X, WS.X, WS.State.Γ, 1, 1, 2, 2, Par, PropB) #compile functions

    aX = @allocated PMFRG.addBL!(
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
    aX += @allocated PMFRG.addBR!(
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
    aXT = @allocated PMFRG.addBLTilde!(WS.Y, WS.X, WS.X, WS.State.Γ, 1, 1, 2, 2, Par, PropB)
    aXT +=
        @allocated PMFRG.addBRTilde!(WS.Y, WS.X, WS.X, WS.State.Γ, 1, 1, 2, 2, Par, PropB)

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
