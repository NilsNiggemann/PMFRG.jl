using PMFRG
test_OneLoopAllocations()
test_DimerFRG(Obsacc = 1e-6)
test_DimerFRG(TwoLoop(),Obsacc = 1e-6,tol = 1e-6) # accuracy of symmetries is finite, given by length of Matsubara sum
test_DimerParquet(tol = 1e-6) # accuracy of symmetries is finite, given by length of Matsubara sum
