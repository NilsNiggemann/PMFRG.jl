"""
Solves FRG as specified for parameters
"""
SolveFRG(Par::PMFRGCore.TwoLoopParams; kwargs...) = launchPMFRG!(
    PMFRGCore.InitializeState(Par),
    PMFRGCore.AllocateSetup(Par),
    PMFRGCore.getDerivTwoLoop!;
    kwargs...,
)
