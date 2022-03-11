using Test,PMFRG
using PMFRG.SpinFRGLattices

"""Exectutes nontrivial symmetry in flow equations for Heisenberg dimer. Local and nonlocal b vertex are equal: Γb_11 = Γb_12!"""
function DimerTest(Method = OneLoop())
    Par = Params(System = getPolymer(2),N=20,T=0.75,accuracy = 1e-3,usesymmetry = true,Lam_min = 0.)
    SolP,ObsPt = SolveFRG(Par,Method)
    vb = SolP.u[end].x[4]
    vb1 =  @view vb[1,:,:,:] 
    vb2 =  @view vb[2,:,:,:] 
    @testset "Dimer: Γb_11 == Γb_12" begin
        for i in eachindex(vb1,vb2)
            @test isapprox(vb1[i],vb2[i],atol = 1e-15)
        end
    end
end