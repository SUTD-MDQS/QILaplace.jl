using QILaplace
using Documenter, Literate

DocMeta.setdocmeta!(QILaplace, :DocTestSetup, :(using QILaplace); recursive=true)

# Generate tutorials
const TUTORIALS_DIR = joinpath(@__DIR__, "src", "tutorials")
const TUTORIALS = [
    "dft",
]

# for tutorial in TUTORIALS
#     file_path = joinpath(TUTORIALS_DIR, "$(tutorial).jl")
#     if isfile(file_path)
#         Literate.markdown(
#             file_path,
#             TUTORIALS_DIR;
#             execute=true,
#         )
#     else
#         @warn "Tutorial source file not found: $file_path"
#     end
# end

makedocs(;
    modules=[
        QILaplace,
        QILaplace.Mps,
        QILaplace.Mpo,
        QILaplace.RSVD,
        QILaplace.ApplyMPO,
        QILaplace.SignalConverters,
        QILaplace.QFTTransform,
        QILaplace.DTTransform,
        QILaplace.ZTTransformer
    ],
    authors="Gauthameshwar S., Noufal Jaseem",
    sitename="QILaplace.jl",
    format=Documenter.HTML(;
        canonical="https://SUTD-MDQS.github.io/QILaplace.jl",
        edit_link="master",
        assets=[
            "assets/sidebar-logo.css",
            "assets/dark-theme.css",
        ],
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        # "Tutorials" => [
        #     "Discrete Fourier Transform" => "tutorials/dft.md",
        #     "Damping Transform" => "tutorials/dt.md",
        #     "Discrete Laplace Transform" => "tutorials/zt.md",
        # ],
        "Benchmarking Runs" => "benchmarking.md",
        "API" => "api.md",
    ],
    checkdocs=:none,
)

deploydocs(;
    repo="github.com/SUTD-MDQS/QILaplace.jl",
    devbranch="master",
)
