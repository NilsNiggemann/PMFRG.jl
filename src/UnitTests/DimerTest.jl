"""Exectutes nontrivial symmetry in flow equations for Heisenberg dimer. Local and nonlocal b vertex are equal: Γb_11 = Γb_12!"""
function test_runFRG(Method = OneLoop(),System = getPolymer(2))
    # Par = Params(getPolymer(2),Method,N=24,Ngamma = 24,T=0.5,accuracy = 1e-3,usesymmetry = false,Lam_min = 0.,MinimalOutput = true,lenIntw = 60)
    Par = BenchmarkingParams(Method,System)
    tempFolder = "temp_PMFRG_test"
    mainFile = joinpath(tempFolder,"temp_main.h5")
    CheckPoints = joinpath(tempFolder,"Checkpoints.h5")
    
    SolP,ObsPt = SolveFRG(Par,MainFile = mainFile,CheckpointDirectory = CheckPoints)
    
    println("cleaning up... deleting ",mainFile, " and ", CheckPoints)
    rm(tempFolder,recursive = true)
    Γa = SolP.u[end].x[3]
    Γb = SolP.u[end].x[4]
    Γc = SolP.u[end].x[5]
    return Γa,Γb,Γc,ObsPt.saveval[end],Par
end


function test_DimerFRG(Method = OneLoop();Obsacc = 1e-14,kwargs...)
    Γa,Γb,Γc,Obs,Par = test_runFRG(Method,getPolymer(2))
    test_FRGResults(Γa,Γb,Γc,Obs,Par;Obsacc = Obsacc,kwargs...)
    
    println("Testing whether local and nonlocal b vertex are equal on dimer Γb_11 = Γb_12")
    test_Gammab_Dimer(Γb;kwargs...)

end

function test_SquagomeFRG(Method = OneLoop();Obsacc = 1e-14,kwargs...)
    SysFunc = SquareKagome.getMirrorSquareKagome
    Γa,Γb,Γc,Obs,Par = test_runFRG(Method,SysFunc(4,1,0.2))
    test_FRGResults(Γa,Γb,Γc,Obs,Par,(Method,SysFunc);Obsacc = Obsacc,kwargs...)
    
end

function test_FRGResults(Γa,Γb,Γc,Obs,Par,Method = PMFRG.getPMFRGMethod(Par);Obsacc=1e-14,kwargs...)
    test_nonzero(Γa,Γb,Γc,Par.System.couplings)

    @testset verbose = true "Testing frequency symmetries" begin
        println("Testing onsite Γa_ii vertex")
        test_Gammaa_onsite(Γa;kwargs...)

        println("Testing t ↔ u symmetries")
        test_tu_symmetries(Γa,Γb,Γc;kwargs...)    
    end

    test_Observables(Method,Obs,Obsacc = Obsacc)

    return
end

function test_runDimerParquet()
    ParquetLambda = 0.
    Par = BenchmarkingParams(Parquet())
    Sol,Obs = SolveParquet(Par,ParquetLambda)
    return Sol,Obs,Par
end

function test_DimerParquet(;kwargs...)
    Sol,Obs,Par = test_runDimerParquet()
    Γa,Γb,Γc = Sol.State.Γ.a,Sol.State.Γ.b,Sol.State.Γ.c
    test_BareBubbles(Sol;kwargs...)
    test_BubbleSymmetries(Sol;kwargs...)
    test_FRGResults(Γa,Γb,Γc,Obs,Par;kwargs...)
end

function test_Observables(Method,Obs;Obsacc=1e-14)
    obs_ex = example_Obs(Method)
    function test(ObsName)
        O1 = getproperty(Obs,ObsName)
        O2 = getproperty(obs_ex,ObsName)
        for i in eachindex(O2)
            @test O1[i] ≈ O2[i] atol = Obsacc
        end
    end

    println("Observables: ", Obs)
    @testset "Testing Susceptibility" begin
        test(:Chi)
    end
    
    @testset "Testing γ" begin
        test(:gamma)
    end
    
    @testset "Testing max(Γ)" begin
        @testset "Γa" begin
            test(:MaxVa)
        end
        @testset "Γb" begin
            test(:MaxVb)
        end
        @testset "Γc" begin
            test(:MaxVc)
        end
    end
end

test_Observables(Method::Parquet,Obs;Obsacc=1e-14) = println("Observables check not implemented for parquet")

function test_nonzero(Γa,Γb,Γc,couplings)
    @testset "non-trivial Vertices" begin
        @testset "Γa" begin
            @test maximum(abs,Γa) >0.
        end
        @testset "Γb" begin
            @test maximum(abs,Γb) > 0.
        end
        @testset "Γc" begin
            @test maximum(abs,Γc) > maximum(abs,couplings)
        end
    end
end

function test_tu_symmetries(Γa,Γb,Γc;kwargs...)
    test_tu_symmetry_ab(Γa,"Γa";kwargs...)
    test_tu_symmetry_ab(Γb,"Γb";kwargs...)
    test_tu_symmetry_c(Γa,Γb,Γc;kwargs...)
end

function test_Gammab_Dimer(Γb::AbstractArray;tol = 1e-14)
    vb1 =  @view Γb[1,:,:,:] 
    vb2 =  @view Γb[2,:,:,:]
    @testset "Dimer: Γb_11 == Γb_12" begin
        Failures = 0
        for i in eachindex(vb1,vb2)
            if !isapprox(vb1[i],vb2[i], atol = tol) 
                Failures += 1
            end
        end
        @test Failures == 0
    end
end

function test_Gammaa_onsite(Γa::AbstractArray,OnsiteBonds=[1];tol = 1e-14)  
    @testset "Γa_ii(stu) == -Γa_ii(uts)" begin
        Failures = 0
        for Rij in OnsiteBonds, s in axes(Γa,2), t in axes(Γa,3), u in axes(Γa,4) 
            if !isapprox(Γa[Rij,s,t,u] , -Γa[Rij,u,t,s], atol = tol) 
                Failures += 1
            end
        end
        @test Failures == 0
    end

    @testset "Γa_ii(stu) == Γa_ii(tus)" begin
        Failures = 0
        for Rij in OnsiteBonds, s in axes(Γa,2), t in axes(Γa,3), u in axes(Γa,4) 
            if !isapprox(Γa[Rij,s,t,u] , +Γa[Rij,t,u,s], atol = tol) 
                Failures += 1
            end
        end
        @test Failures == 0
    end
end

function test_tu_symmetry_ab(Γ::AbstractArray,Name;tol = 1e-14) 
    
    @testset "$(Name)_ij(stu) == -$(Name)_ij(sut)" begin
        Failures = 0
        for Rij in axes(Γ,1), s in axes(Γ,2), t in axes(Γ,3), u in axes(Γ,4) 
            if !isapprox(Γ[Rij,s,t,u] , -Γ[Rij,s,u,t], atol = tol) 
                Failures += 1
            end
        end
        @test Failures == 0

    end
end

function test_tu_symmetry_c(Γa::AbstractArray,Γb::AbstractArray,Γc::AbstractArray,V = "Γ";tol = 1e-14)
    Failures = 0
    @testset "$(V)c_ij(stu) == (-$(V)a_ij + $(V)b_ij + $(V)c_ij)(sut)" begin
        for Rij in axes(Γc,1), s in axes(Γc,2), t in axes(Γc,3), u in axes(Γc,4) 
            if !isapprox(Γc[Rij,s,t,u] , (-Γa[Rij,s,u,t] +Γb[Rij,s,u,t] +Γc[Rij,s,u,t]), atol = tol) 
                Failures += 1
            end
        end
        @test Failures == 0
    end

end

