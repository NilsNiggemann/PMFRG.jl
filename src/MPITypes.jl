struct UseMPI end

struct MPIOneLoopParams{F<:AbstractFloat,G<:Geometry} <: PMFRG.AbstractOneLoopParams
    System::G
    NumericalParams::NumericalParams{F}
    Options::OptionParams
end

Params(System::PMFRG.Geometry, ::UseMPI; kwargs...) =
    MPIOneLoopParams(System, NumericalParams(; kwargs...), OptionParams(; kwargs...)) # Todo: make this error when an unknown kwarg is given!
