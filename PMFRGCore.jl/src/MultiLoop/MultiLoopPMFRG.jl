"""
includes MultiLoop corrections to PMFRG.
"""

include("MultiLoopTypes.jl")
include("Parquet.jl")

generateFileName(Par::MultiLoopParams, arg::String = "") =
    _generateFileName(Par, "_l$(Par.l)" * arg)
generateFileName(Par::ParquetParams, arg::String = "") = _generateFileName(Par, "_p" * arg)
