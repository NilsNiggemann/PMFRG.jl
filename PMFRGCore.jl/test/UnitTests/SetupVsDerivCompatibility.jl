function test_setup_deriv_compatibility()

    @testset "getDeriv! can be called with the state and the setup allocated" begin
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

        State = PMFRGCore.InitializeState(par)
        Deriv = similar(State)

        setup = PMFRGCore.AllocateSetup(par)
        Lam = 1.0

        @testset "Plugging all together..." begin
            @test isa(PMFRGCore.getDeriv!(Deriv, State, setup, Lam), Any)
        end
    end
end
