using PMFRG
using Test
using Serialization
thisdir = dirname(@__FILE__)
data = deserialize(joinpath(thisdir,"PMFRG.getXBubble.data"))

include("PMFRG.getXBubble.common.jl")

