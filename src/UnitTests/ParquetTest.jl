function test_SDE(Par::ParquetParams)
    Workspace = SetupParquet(Par)
    Lam = 0.
    @unpack OldState,State,Γ0,X,B0,BX,Par,Buffer = Workspace

    @inline Prop(x,nw) = 1/6*iG_(State.γ,x,Lam,nw,Par.NumericalParams.T)

    getProp! = constructPropagatorFunction(Workspace,Lam)

    gamma1 = copy(State.γ)
    gamma2 = copy(State.γ)

    computeLeft2PartBubble!(B0,Γ0,Γ0,State.Γ,getProp!,Par,Buffer)

    compute1PartBubble!(gamma1,B0,Prop,Par)
    compute1PartBubble_BS!(gamma2,State.Γ,Γ0,Prop,Par)
    @testset "SDE gamma" begin
        for i in eachindex(gamma1,gamma2)
            @test gamma1[i] ≈ gamma2[i]
        end
    end
    State.Γ.c |> minimum
end

function test_FirstBSEIteration!(Workspace,Lam)
    @unpack OldState,State,Γ0,X,B0,BX,Par,Buffer = Workspace

    @inline Prop(x,nw) = 1/6*iG_(State.γ,x,Lam,nw,Par.NumericalParams.T)

    getProp! = constructPropagatorFunction(Workspace,Lam)

    computeLeft2PartBubble!(B0,Γ0,Γ0,State.Γ,getProp!,Par,Buffer)
end

function test_SDE_FP(Par::ParquetParams)
    Workspace1 = SetupParquet(Par)
    Workspace2 = SetupParquet(Par)
    Lam = 0.
    test_FirstBSEIteration!(Workspace1,Lam)
    test_FirstBSEIteration!(Workspace2,Lam)

    iterateSDE!(Workspace1,Lam)
    iterateSDE_FP!(Workspace2,Lam)
    @testset "SDE FixedPoint" begin
        State1 = Workspace1.State
        State2 = Workspace2.State
        for i in eachindex(State1.γ,State2.γ)
            @test State1.γ[i] ≈ State2.γ[i]
        end
        for field in (:a,:b,:c)
            V1 = getproperty(State1.Γ,field)
            V2 = getproperty(State2.Γ,field)
            for i in eachindex(V1,V2)
                @test V1[i] ≈ V2[i]
            end
        end

    end
end