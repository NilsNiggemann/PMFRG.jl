using PMFRG, PMFRG.UnitTests
test_OneLoopAllocations()
test_DimerFRG()
test_DimerFRG(TwoLoop(),tol = 1e-5) # accuracy of symmetries is finite, given by length of Matsubara sum
