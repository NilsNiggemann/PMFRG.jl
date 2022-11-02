using SpinFRGLattices, PMFRG
using SpinFRGLattices.SquareLattice
##
NLen = 5
J1 = 1
J2 = 0.1
couplings = [J1,J2] 
System = getSquareLattice(NLen,couplings) 

Trange = 0.3:0.05:1.5
for T in Trange
    Par = Params(
        System,
        OneLoop(),
        T = T,
        N = 10,
        accuracy = 1e-8,
        MinimalOutput = true,
    )

    mainFile = "temp/"*PMFRG.generateFileName(Par,"_testFile_TSweep")

    Solution,saved_values = SolveFRG(Par,MainFile = mainFile,method = DP5())
end