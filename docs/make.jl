using QILaplace
using Documenter

DocMeta.setdocmeta!(QILaplace, :DocTestSetup, :(using QILaplace); recursive=true)

makedocs(;
    modules=[QILaplace],
    authors="Gauthameshwar S., Noufal Jaseem",
    sitename="QILaplace.jl",
    format=Documenter.HTML(;
        canonical="https://SUTD-MDQS.github.io/QILaplace.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/SUTD-MDQS/QILaplace.jl",
    devbranch="master",
)
