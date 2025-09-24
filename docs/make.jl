using SpatialMultitaper, Documenter, DocumenterCitations, DocumenterMermaid, Literate

## Literate preprocessing, maybe move to a separate script later for faster builds

# to run with LiveServer and avoid infinite loops, use
# servedocs(literate_dir=joinpath("docs","literate","tutorials"),skip_dir = joinpath("docs","src","tutorials"))
# adapting `tutorials` to whatever subdir you are working on

LITERATE_INPUT = joinpath(@__DIR__, "literate")
LITERATE_OUTPUT = joinpath(@__DIR__, "src")

for dir_path in filter(isdir, readdir(joinpath(@__DIR__, "literate"), join = true))
    dirname = basename(dir_path)

    for (root, _, files) in walkdir(dir_path), file in files
        if occursin("dump", root)
            @warn "Skipping dump file $file"
            continue
        end
        # ignore non julia files
        splitext(file)[2] == ".jl" || continue
        # full path to a literate script
        ipath = joinpath(root, file)
        # generated output path
        opath = splitdir(replace(ipath, LITERATE_INPUT => LITERATE_OUTPUT))[1]
        # generate the markdown file calling Literate
        Literate.markdown(ipath, opath)
    end
end

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "citation.bib");
    style = :numeric
)

DocMeta.setdocmeta!(
    SpatialMultitaper, :DocTestSetup, :(using SpatialMultitaper); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [file
                        for file in readdir(joinpath(@__DIR__, "src"))
                        if
                        file != "index.md" && splitext(file)[2] == ".md"]

makedocs(;
    modules = [SpatialMultitaper],
    authors = "Jake P Grainger",
    repo = "https://github.com/JakeGrainger/SpatialMultitaper.jl/blob/{commit}{path}#{line}",
    sitename = "SpatialMultitaper.jl",
    format = Documenter.HTML(;
        canonical = "https://JakeGrainger.github.io/SpatialMultitaper.jl"),
    plugins = [bib],
    pages = ["index.md";
             "tutorials" => ["basics" => "tutorials/basic_estimation.md",
                 "tapers" => "tutorials/tapers.md"]
             numbered_pages]
)

deploydocs(; repo = "github.com/JakeGrainger/SpatialMultitaper.jl")
