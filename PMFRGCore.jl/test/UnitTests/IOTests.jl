function test_saving(
    FileDirectory,
    Method = OneLoop(),
    GeometryGenerator = SquareKagome.getSquareKagome,
)
    FileDirectory = PMFRG.UniqueDirName(FileDirectory)
    mkpath(FileDirectory)
    Par = BenchmarkingParams(Method, GeometryGenerator(4))
    State = PMFRG.InitializeState(Par)
    Lam = Par.NumericalParams.Lam_min + 0.01
    saved_values = DiffEqCallBacks.SavedValues(Float64, PMFRG.Observables)

    obs1 = PMFRG.getObservables(PMFRG.Observables, State, Par.NumericalParams.Lam_max, Par)
    obs2 = PMFRG.getObservables(PMFRG.Observables, State, Lam, Par)

    push!(saved_values.t, Par.NumericalParams.Lam_max)
    push!(saved_values.saveval, obs1)

    push!(saved_values.t, Lam)
    push!(saved_values.saveval, obs2)

    Filename = PMFRG.saveCurrentState(FileDirectory, State, saved_values, Lam, Par)
    return Filename
end

function test_loading(
    Filename::String,
    GeometryGenerator::Function = SquareKagome.getSquareKagome;
    kwargs...,
)
    Obs = PMFRG.readObservables(Filename)
    Sol, saved_values = SolveFRG_Checkpoint(
        Filename,
        readGeometry(Filename, GeometryGenerator);
        Params = (MinimalOutput = true,),
        kwargs...,
    )
    @testset "Observables extended" begin
        @test issubset(Obs.t, saved_values.t)
    end
end

function test_IO(Method = OneLoop(), GeometryGenerator = SquareKagome.getSquareKagome)
    @testset "FileIO" begin
        tempFolder = joinpath("temp_PMFRG_test")
        CheckpointDir = joinpath("temp_PMFRG_test", "Checkpoints")
        Filename = test_saving(CheckpointDir, Method, GeometryGenerator)

        Par = PMFRG.readParams(Filename, GeometryGenerator)
        @testset "Testing loopOrder from Params" begin
            @test PMFRG.getPMFRGMethod(Par) == Method
        end
        test_loading(Filename, GeometryGenerator)
        rm(tempFolder, recursive = true)
        println("cleaned up... deleted ", tempFolder)
    end
end
