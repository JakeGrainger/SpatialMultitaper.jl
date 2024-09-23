using SpatialMultitaper
using Documenter

DocMeta.setdocmeta!(SpatialMultitaper, :DocTestSetup, :(using SpatialMultitaper); recursive=true)

makedocs(;
    modules=[SpatialMultitaper],
    authors="Jake Grainger",
    repo="https://github.com/JakeGrainger/SpatialMultitaper.jl/blob/{commit}{path}#{line}",
    sitename="SpatialMultitaper.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JakeGrainger.github.io/SpatialMultitaper.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JakeGrainger/SpatialMultitaper.jl",
    devbranch="main",
)
