"""
includes MultiLoop corrections to PMFRG.
"""

include("MultiLoopTypes.jl")
include("Flowequations.jl")
include("Parquet.jl")

export SolveParquet
generateFileName(Par::MultiLoopParams,arg::String = "") = generateFileName(Par,string("_",Method.l,arg))
