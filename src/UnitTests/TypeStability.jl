
"""gets empty Workspace."""
function getEmptyWorkspace(Par::PMFRG.OneLoopParams)
    State = PMFRG.InitializeState(Par)
    Deriv = similar(State)
    setup = PMFRG.AllocateSetup(Par)
    (X,Buffs,Par) = setup
    PMFRG.OneLoopWorkspace(Deriv,State,X,Buffs,Par)
end

"""Checks whether vertex derivative has type instabilities"""
function test_OneLoopAllocations(Par=Params(getPolymer(2)))
    WS = getEmptyWorkspace(Par)
    # eval("@code_warntype PMFRG.getXBubble!(WS,0.)")
    # @code_warntype PMFRG.getXBubble!(WS,0.)
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

    Npairs =Par.System.Npairs
    N = Par.NumericalParams.N

    WS.X.a .= rand(1:0.1:10, Npairs,N,N,N)
    WS.X.b .= rand(1:0.1:10, Npairs,N,N,N)
    WS.X.c .= rand(1:0.1:10, Npairs,N,N,N)
    WS.X.Ta .= rand(1:0.1:10, Npairs,N,N,N)
    WS.X.Ta .= rand(1:0.1:10, Npairs,N,N,N)
    WS.X.Ta .= rand(1:0.1:10, Npairs,N,N,N)

    WS.State.Γ.a .= rand(1:0.1:10, Npairs,N,N,N)
    WS.State.Γ.b .= rand(1:0.1:10, Npairs,N,N,N)
    WS.State.Γ.c .= rand(1:0.1:10, Npairs,N,N,N)
    PropB .= 1.

    aX = @allocated PMFRG.addBL!(WS.Y,WS.X,WS.X,WS.State.Γ,1,1,2,2,Par,PropB,VB,XB)
    aX += @allocated PMFRG.addBR!(WS.Y,WS.X,WS.X,WS.State.Γ,1,1,2,2,Par,PropB,VB,XB)
    aXT = @allocated PMFRG.addBLTilde!(WS.Y,WS.X,WS.X,WS.State.Γ,1,1,2,2,Par,PropB)
    aXT += @allocated PMFRG.addBRTilde!(WS.Y,WS.X,WS.X,WS.State.Γ,1,1,2,2,Par,PropB)

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
