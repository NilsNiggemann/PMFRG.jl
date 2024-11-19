using Test

using MPI, PencilArrays

using PMFRGCore
using SpinFRGLattices


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



function smoketest_function_compatibilities_MPI()
    @testset verbose = true "MPI-related Smoke tests " begin
        @testset "Checking that PMFRGCoreMPIExt is loaded" begin
            @test !isnothing(Base.get_extension(PMFRGCore, :PMFRGCoreMPIExt))
        end
        @testset verbose = true "compat between funcs using/creating State/setup" begin
            test_state_init(par)
            test_state_setup_deriv_mpi(par)
            test_state_getObservables_mpi(par)
        end
    end


end

function test_state_init(par)
    State = PMFRGCore.InitializeState(par, PMFRGCore.UseMPI())
    @test isa(State, PencilArray)
    (; couplings) = par.System

    @test sum(State) == -sum(couplings) * par.NumericalParams.N^3

end


function test_state_setup_deriv_mpi(par)
    @testset "getDeriv! can be called with the state and the setup allocated" begin
        State = PMFRGCore.InitializeState(par, PMFRGCore.UseMPI())
        Deriv = similar(State)

        setup = PMFRGCore.AllocateSetup(par, PMFRGCore.UseMPI())
        Lam = 1.0

        @testset "Plugging all together..." begin
            @test isa(PMFRGCore.getDeriv!(Deriv, State, setup, Lam), Any)
        end
    end
end


function test_state_getObservables_mpi(par)
    @testset "getObservables can be called with the state allocated" begin
        State = PMFRGCore.InitializeState(par, PMFRGCore.UseMPI())

        @test isa(PMFRGCore.getObservables(PMFRGCore.Observables, State, 0.5, par), Any)
    end
end


MPI.Init()
@testset verbose = true "Tests that need MPI to run" begin
    @testset "type of ParallelizationScheme in setup with UseMPI is UseMPI" begin
        setup = PMFRGCore.AllocateSetup(par, PMFRGCore.UseMPI())
        @test isa(setup.ParallelizationScheme, PMFRGCore.UseMPI)
    end
    smoketest_function_compatibilities_MPI()
end
MPI.Finalize()
