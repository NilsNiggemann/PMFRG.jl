"""
Solves FRG as specified for parameters
"""
SolveFRG(Par::TwoLoopParams; kwargs...) =
    launchPMFRG!(InitializeState(Par), AllocateSetup(Par), getDerivTwoLoop!; kwargs...)
