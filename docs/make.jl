using Pkg
Pkg.activate()
using Documenter
Pkg.activate("../")
push!(LOAD_PATH,"../src/")

using PMFRG
makedocs(sitename = "PMFRG.jl",modules  = [PMFRG],pages=["Home" => "index.md"])
deploydocs(;repo="git@gitlabph.physik.fu-berlin.de:niggeni/PMFRG.jl",)
# deploydocs(;repo="github.com/NilsNiggemann/PMFRG.jl",)