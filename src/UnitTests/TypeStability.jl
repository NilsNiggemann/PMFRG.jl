
"""gets empty Workspace."""
function getEmptyWorkspace(Par::PMFRG.OneLoopParams)
    State = PMFRG.InitializeState(Par)
    Deriv = similar(State)
    setup = PMFRG.AllocateSetup(Par)
    (X,Buffs,Par) = setup
    PMFRG.OneLoopWorkspace(Deriv,State,X,Buffs,Par)
end

convertToSMatrix(M::MMatrix{N,N,T,NN}) where {N,T,NN} = SMatrix{N,N,T,NN}(M)

"""Checks whether vertex derivative has type instabilities"""
function test_OneLoopAllocations(Par=Params(getPolymer(2)))
    WS = getEmptyWorkspace(Par)
    # eval("@code_warntype PMFRG.getXBubble!(WS,0.)")
    # @code_warntype PMFRG.getXBubble!(WS,0.)

    SProps = convertToSMatrix(WS.Buffer.Props[1])
    VBuff = WS.Buffer.Vertex[1]
    test_OneLoopAllocations(Par,WS,SProps,VBuff)
end
function test_OneLoopAllocations(Par,WS,SProps,VBuff)
    addX!(WS,1,1,2,2,SProps,VBuff) # compile functions
    addXTilde!(WS,1,1,2,2,SProps) # compile functions

    aX = @allocated addX!(WS,1,1,2,2,SProps,VBuff)
    aXT = @allocated addXTilde!(WS,1,1,2,2,SProps)
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
    
    PropB =convertToSMatrix(WS.Buffer.Props[1])

    VB = WS.Buffer.Vertex[1]
    XB =  WS.Buffer.X[1]

    test_TwoLoopAllocations(Par,WS,PropB,VB,XB)
end
function test_TwoLoopAllocations(Par,WS,PropB,VB,XB)
    addBL!(WS.Y,WS.X,WS.X,WS.State.Γ,1,1,2,2,Par,PropB,VB,XB) #compile functions
    addBR!(WS.Y,WS.X,WS.X,WS.State.Γ,1,1,2,2,Par,PropB,VB,XB) #compile functions
    addBLTilde!(WS.Y,WS.X,WS.X,WS.State.Γ,1,1,2,2,Par,PropB) #compile functions
    addBRTilde!(WS.Y,WS.X,WS.X,WS.State.Γ,1,1,2,2,Par,PropB) #compile functions
    
    aX = @allocated addBL!(WS.Y,WS.X,WS.X,WS.State.Γ,1,1,2,2,Par,PropB,VB,XB)
    aX += @allocated addBR!(WS.Y,WS.X,WS.X,WS.State.Γ,1,1,2,2,Par,PropB,VB,XB)
    aXT = @allocated addBLTilde!(WS.Y,WS.X,WS.X,WS.State.Γ,1,1,2,2,Par,PropB)
    aXT += @allocated addBRTilde!(WS.Y,WS.X,WS.X,WS.State.Γ,1,1,2,2,Par,PropB)
    
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
