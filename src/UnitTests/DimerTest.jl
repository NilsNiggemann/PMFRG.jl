"""Exectutes nontrivial symmetry in flow equations for Heisenberg dimer. Local and nonlocal b vertex are equal: Γb_11 = Γb_12!"""
function test_DimerFRG(Method = OneLoop();Obsacc = 1e-6,kwargs...)
    Par = Params(getPolymer(2),Method,N=24,Ngamma = 24,T=0.5,accuracy = 1e-3,usesymmetry = false,Lam_min = 0.,MinimalOutput = true,lenIntw = 60)

    tempFolder = "temp_PMFRG_test2"
    mainFile = joinpath(tempFolder,"temp_main.h5")
    CheckPoints = joinpath(tempFolder,"Checkpoints.h5")
    
    SolP,ObsPt = SolveFRG(Par,MainFile = mainFile,CheckpointDirectory = CheckPoints)

    println("cleaning up... deleting ",mainFile, " and ", CheckPoints)
    rm(tempFolder,recursive = true)

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

    test_Observables(Method,ObsPt.saveval[end],Obsacc = Obsacc)

    return
end

function test_Observables(Method,Obs;Obsacc=1e-6)
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

function example_Obs(Method::OneLoop)
    PMFRG.Observables([0.409961674309049, -0.200041163060498], [0.367090761453024 0.128045613399846 0.086744514258760 0.060408397243958 0.046089629935719 0.037224838436830 0.031123238584940 0.026623695820601 0.023146972150518 0.020378546539016 0.018118952133753 0.016242668556776 0.014654076476321 0.013290063763770 0.012095764689258 0.011036461682002 0.010077688811575 0.009199596223172 0.008378589894997 0.007602730223307 0.006854346043725 0.006127144022231 0.005429363295683 0.001190178969660], [-0.095054800124967], [0.035369300845790, 1.789921595488008], [0.755459278321820, 0.755459278321820], [0.755459278321820, 2.689104516815582])
end

function example_Obs(Method::TwoLoop)
    PMFRG.Observables([0.384798093987635, -0.180469046203119], [0.482956311303588 0.222352984832864 0.154573855716267 0.109897378628719 0.084917553190298 0.069199964487509 0.058308013987969 0.050266564388793 0.044063205870330 0.039131077027183 0.035112384954753 0.031778170059306 0.028961474004695 0.026549342793794 0.024450131950865 0.022602237011145 0.020950360390247 0.019459109883070 0.018092609825221 0.016830289679813 0.015649404101204 0.014548321178442 0.013659541068861 -0.005699647497969], [-0.100154313528458], [0.222988595275207, 2.111032725409252], [1.106789531343003, 1.106789531371678], [1.106789531343003, 2.660858595690006])
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

