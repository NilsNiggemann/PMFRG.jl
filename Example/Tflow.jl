using SpinFRGLattices, PMFRG
Par = Params(getPolymer(2),N=10,T_min = 0.4,T_max =exp(7),accuracy = 1e-5)
SolveFRG(Par)