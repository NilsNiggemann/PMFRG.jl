using Test
using PMFRG
using SpinFRGLattices.SquareLattice


function _getsystem()
    NLen = 5
    J1 = 1
    J2 = 0.1
    couplings = [J1, J2]
    getSquareLattice(NLen, couplings)
end


function test_params_wrong_kwargs_oneloop()
    @test_throws ArgumentError Params(
        _getsystem(),
        OneLoop(),
        T = 0.5,
        N = 10,
        accuracy = 1e-3,
        someRandomKwarg = "impostor",
    )

    @test_throws "someRandomKwarg" Params(
        _getsystem(),
        OneLoop(),
        T = 0.5,
        N = 10,
        accuracy = 1e-3,
        someRandomKwarg = "impostor",
    )


end

function test_params_wrong_kwargs_twoloop()
    @test_throws ArgumentError Params(
        _getsystem(),
        TwoLoop(),
        T = 0.5,
        N = 10,
        accuracy = 1e-3,
        someRandomKwarg = "impostor",
    )

    @test_throws "someRandomKwarg" Params(
        _getsystem(),
        TwoLoop(),
        T = 0.5,
        N = 10,
        accuracy = 1e-3,
        someRandomKwarg = "impostor",
    )

end

function test_params_wrong_kwargs_multiloop()
    @test_throws ArgumentError Params(
        _getsystem(),
        MultiLoop(3),
        T = 0.5,
        N = 10,
        accuracy = 1e-3,
        someRandomKwarg = "impostor",
    )

    @test_throws "someRandomKwarg" Params(
        _getsystem(),
        MultiLoop(3),
        T = 0.5,
        N = 10,
        accuracy = 1e-3,
        someRandomKwarg = "impostor",
    )

end

function test_params_wrong_kwargs_parquet()
    @test_throws ArgumentError Params(
        _getsystem(),
        Parquet(),
        T = 0.5,
        someRandomKwarg = "impostor",
    )

    @test_throws "someRandomKwarg" Params(
        _getsystem(),
        Parquet(),
        T = 0.5,
        someRandomKwarg = "impostor",
    )
end


function test_params_wrong_kwargs()
    @testset verbose = true begin
        test_params_wrong_kwargs_oneloop()
        test_params_wrong_kwargs_twoloop()
        test_params_wrong_kwargs_multiloop()
        test_params_wrong_kwargs_parquet()
    end
    return
end
