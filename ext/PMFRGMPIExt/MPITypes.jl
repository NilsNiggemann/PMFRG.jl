
struct UseMPI end

struct MPIOneLoopParams{F<:AbstractFloat,G<:PMFRG.Geometry} <: PMFRG.OneLoopParams
    System::G
    NumericalParams::PMFRG.NumericalParams{F}
    Options::PMFRG.OptionParams
end

PMFRG.Params(System::PMFRG.Geometry, ::UseMPI; kwargs...) = MPIOneLoopParams(
    System,
    PMFRG.NumericalParams(; kwargs...),
    PMFRG.OptionParams(; kwargs...),
) # Todo: make this error when an unknown kwarg is given!
