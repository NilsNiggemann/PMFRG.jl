function Testrun()
    System = getPolymer(2)
    Par = Params(System = System,N=10,T=1.5,MinimalOutput=true,usesymmetry = true,accuracy= 1E-5,Ngamma = 100)
    println("precompiling: ")
    Solution,saved_values = SolveFRG(Par,method = DP5())
    # Solution,saved_values = SolveFRG(Par,method = Vern7())
    return Par,Solution,saved_values
end
