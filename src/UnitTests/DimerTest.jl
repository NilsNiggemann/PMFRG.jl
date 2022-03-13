"""Exectutes nontrivial symmetry in flow equations for Heisenberg dimer. Local and nonlocal b vertex are equal: Γb_11 = Γb_12!"""
function test_DimerFRG(Method = OneLoop();kwargs...)
    Par = Params(getPolymer(2),Method,N=24,T=0.5,accuracy = 1e-3,usesymmetry = false,Lam_min = 0.,MinimalOutput = true,lenIntw = 60)

    tempFolder = "temp_PMFRG_test"
    mainFile = joinpath(tempFolder,"temp_main.h5")
    CheckPoints = joinpath(tempFolder,"Checkpoints.h5")
    
    SolP,ObsPt = SolveFRG(Par,MainFile = mainFile,CheckpointDirectory = CheckPoints)

    println("cleaning up... deleting ",mainFile, " and ", CheckPoints)

    Chi = ObsPt.saveval[end].Chi
    @testset "Testing Susceptibility" begin
        println("χ = ", Chi)
        @test Chi[1] >0. &&Chi[1]<1.
        @test Chi[2] <0. &&Chi[1]>-1.
    end
    
    Γa = SolP.u[end].x[3]
    Γb = SolP.u[end].x[4]
    Γc = SolP.u[end].x[5]
    test_nonzero(Γa,Γb,Γc,Par.System.couplings)

    @testset verbose = true "Testing frequency symmetries" begin
        println("Testing onsite Γa_ii vertex")
        test_Gammaa_onsite(Γa;kwargs...)

        println("Testing t ↔ u symmetries")
        test_tu_symmetries(Γa,Γb,Γc;kwargs...)    
    end

    println("Testing whether local and nonlocal b vertex are equal on dimer Γb_11 = Γb_12")
    test_Gammab_Dimer(Γb;kwargs...)


    rm(tempFolder,recursive = true)
    return
end

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

function test_Gammab_Dimer(Γb::AbstractArray;tol = 1e-15)
    vb1 =  @view Γb[1,:,:,:] 
    vb2 =  @view Γb[2,:,:,:] 
    @testset "Dimer: Γb_11 == Γb_12" begin
        for i in eachindex(vb1,vb2)
            @test vb1[i] ≈ vb2[i] atol = tol
        end
    end
end

function test_Gammaa_onsite(Γa::AbstractArray,OnsiteBonds=[1];tol = 1e-15)  
    @testset "Γa_ii(stu) == -Γa_ii(uts)" begin
        for Rij in OnsiteBonds, s in axes(Γa,2), t in axes(Γa,3), u in axes(Γa,4) 
            va = +Γa[Rij,s,t,u]
            @test va ≈ -Γa[Rij,u,t,s] atol = tol
        end
    end
    @testset "Γa_ii(stu) == Γa_ii(tus)" begin
        for Rij in OnsiteBonds, s in axes(Γa,2), t in axes(Γa,3), u in axes(Γa,4) 
            va = +Γa[Rij,s,t,u]
            @test va ≈ +Γa[Rij,t,u,s] atol = tol
        end
    end
end

function test_tu_symmetry_ab(Γ::AbstractArray,Name;tol = 1e-15)  
    @testset "$(Name)_ij(stu) == -$(Name)_ij(sut)" begin
        for Rij in axes(Γ,1), s in axes(Γ,2), t in axes(Γ,3), u in axes(Γ,4) 
            @test Γ[Rij,s,t,u] ≈ -Γ[Rij,s,u,t] atol = tol
        end
    end
end
function test_tu_symmetry_c(Γa::AbstractArray,Γb::AbstractArray,Γc::AbstractArray;tol = 1e-1)  
    @testset "Γc_ij(stu) == (-Γa_ij + Γb_ij + Γc_ij)(sut)" begin
        
        for Rij in axes(Γc,1), s in axes(Γc,2), t in axes(Γc,3), u in axes(Γc,4) 
        @test Γc[Rij,s,t,u] ≈ (-Γa[Rij,s,u,t] +Γb[Rij,s,u,t] +Γc[Rij,s,u,t]) atol = tol
        end
    end

end

