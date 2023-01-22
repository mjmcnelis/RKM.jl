push!(LOAD_PATH,"../src/")
using Documenter, Revise, RKM
include("pages.jl")
DocMeta.setdocmeta!(RKM, :DocTestSetup, :(using RKM); recursive=true)

makedocs(
    sitename = "RKM.jl",
    modules  = Module[RKM],
    repo     = "https://github.com/mjmcnelis/RKM.jl",
    source   = "src",
    format   = Documenter.HTML(prettyurls = false),
    pages    = pages,
    # strict   = true # See https://juliadocs.github.io/Documenter.jl/stable/lib/public/#Documenter.makedocs
)

deploydocs(
    repo         = "https://github.com/mjmcnelis/RKM.jl.git";
    push_preview = true
)