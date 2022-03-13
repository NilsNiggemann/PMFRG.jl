
"""gets empty Workspace."""
function getEmptyWorkspace(Par::PMFRG.OneLoopParams)
    State = PMFRG.InitializeState(Par)
    Deriv = similar(State)
    setup = PMFRG.AllocateSetup(Par)
    (X,Buffs,Par) = setup
    PMFRG.OneLoopWorkspace(Deriv,State,X,Buffs,Par)
end

"""Checks wheter vertex derivative has type instabilities"""
function test_OneLoopAllocations(Par=Params(getPolymer(2)))
    WS = getEmptyWorkspace(Par)
    # eval("@code_warntype PMFRG.getVertexDeriv!(WS,0.)")
    # @code_warntype PMFRG.getVertexDeriv!(WS,0.)
    aX = @allocated PMFRG.addX!(WS,1,1,2,2,WS.Buffer.Props[1],WS.Buffer.Vertex[1])
    aXT = @allocated PMFRG.addXTilde!(WS,1,1,2,2,WS.Buffer.Props[1])
    test_TypeStability(Par,WS,aX,aXT)
end

function test_TypeStability(Par,WS,aX,aXT)
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
    (X,Y,Buffs,Par) = setup
    PMFRG.TwoLoopWorkspace(Deriv,State,X,Y,Buffs,Par)
end
function test_TwoLoopAllocations(Par=Params(getPolymer(2),TwoLoop()))
    WS = getEmptyWorkspace(Par)

    PropB = WS.Buffer.Props[1]
    VB = WS.Buffer.Vertex[1]
    XB =  WS.Buffer.X[1]

    aX = @allocated PMFRG.addY!(WS,1,1,2,2,PropB,VB,XB)
    aXT = @allocated PMFRG.addYTilde!(WS,1,1,2,2,PropB)

    test_TypeStability(Par,WS,aX,aXT)
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
