using SpatialMultitaper
using Documenter
using DocumenterCitations

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "citation.bib");
    style=:numeric
)

DocMeta.setdocmeta!(SpatialMultitaper, :DocTestSetup, :(using SpatialMultitaper); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
    file for file in readdir(joinpath(@__DIR__, "src")) if
    file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
    modules = [SpatialMultitaper],
    authors = "Jake P Grainger",
    repo = "https://github.com/JakeGrainger/SpatialMultitaper.jl/blob/{commit}{path}#{line}",
    sitename = "SpatialMultitaper.jl",
    format = Documenter.HTML(; canonical = "https://JakeGrainger.github.io/SpatialMultitaper.jl"),
    plugins = [bib],
    pages = ["index.md"; numbered_pages],
)

deploydocs(; repo = "github.com/JakeGrainger/SpatialMultitaper.jl")
