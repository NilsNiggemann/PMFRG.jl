"""
includes MultiLoop corrections to PMFRG.
"""

include("MultiLoopTypes.jl")
include("Parquet.jl")

export SolveParquet
generateFileName(Par::MultiLoopParams,arg::String = "") = _generateFileName(Par,string("_",Method.l,arg))
