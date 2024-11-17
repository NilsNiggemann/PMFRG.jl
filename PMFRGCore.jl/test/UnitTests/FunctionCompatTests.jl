function smoketest_function_compatibilities()
    par = Params(
        getPolymer(2),
        OneLoop(),
        T = 0.5,
        N = 10,
        Ngamma = 10,
        accuracy = 1e-3,
        Lam_min = exp(-30),
        Lam_max = 100.0,
        usesymmetry = false,
        MinimalOutput = true,
        lenIntw = 60,
        lenIntw_acc = 60,
    )


    @testset verbose = true "compat between funcs using/creating State/setup" begin
        test_state_setup_deriv(par)
        test_state_getObservables(par)
    end
end

function test_state_setup_deriv(par)
    @testset "getDeriv! can be called with the state and the setup allocated" begin
        State = PMFRGCore.InitializeState(par)
        Deriv = similar(State)

        setup = PMFRGCore.AllocateSetup(par)
        Lam = 1.0

        @testset "Plugging all together..." begin
            @test isa(PMFRGCore.getDeriv!(Deriv, State, setup, Lam), Any)
        end
    end
end

function test_state_getObservables(par)
    @testset "getObservables can work with the state allocated" begin
        State = PMFRGCore.InitializeState(par)

        @testset "Pluggin all together..." begin
            @test isa(
                PMFRGCore.getObservables(
                    PMFRGCore.Observables,
                    State,
                    PMFRGCore.t_to_Lam(0.5),
                    par,
                ),
                Any,
            )

        end
    end
end
