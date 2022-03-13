
"""gets empty Workspace."""
function getEmptyOneLoopWorkspace(Par)
    State = PMFRG.InitializeState(Par)
    Deriv = similar(State)
    setup = PMFRG.AllocateSetup(Par)
    (X,Buffs,Par) = setup
    PMFRG.OneLoopWorkspace(Deriv,State,X,Buffs,Par)
end

"""Checks wheter vertex derivative has type instabilities"""
function test_Vertex_Derivative(Par=Params(getPolymer(2)))
    WS = getEmptyOneLoopWorkspace(Par)
    # eval("@code_warntype PMFRG.getVertexDeriv!(WS,0.)")
    # @code_warntype PMFRG.getVertexDeriv!(WS,0.)
    aX = @allocated PMFRG.addX!(WS,1,1,2,2,WS.Buffer.Props[1],WS.Buffer.Vertex[1])
    aXT = @allocated PMFRG.addXTilde!(WS,1,1,2,2,WS.Buffer.Props[1])

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
        @testset "zero allocations of addX! " begin
            @test aXT == 0
        end
    end
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
