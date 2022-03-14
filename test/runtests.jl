using PMFRG, PMFRG.UnitTests
test_OneLoopAllocations()
test_DimerFRG(Chiacc = 0.01)
test_DimerFRG(TwoLoop(),Chiacc = 0.01,tol = 1e-5) # accuracy of symmetries is finite, given by length of Matsubara sum
