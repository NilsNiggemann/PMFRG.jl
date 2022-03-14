using PMFRG, PMFRG.UnitTests
test_OneLoopAllocations()
test_DimerFRG(Obsacc = 1e-6)
test_DimerFRG(TwoLoop(),Obsacc = 1e-6,tol = 1e-5) # accuracy of symmetries is finite, given by length of Matsubara sum
