using Documenter:
    DocMeta,
    HTML,
    MathJax3,
    asset,
    deploydocs,
    makedocs
using PlutoStaticHTML

const NOTEBOOK_DIRS = [
                       joinpath(dirname(@__DIR__), "notebooks", "tropical"), 
                       joinpath(dirname(@__DIR__), "notebooks", "yao"), 
                      ]

"""
    build()
Run all Pluto notebooks (".jl" files) in `NOTEBOOK_DIRS`.
"""
function build()
    println("Building notebooks")
    hopts = HTMLOptions(; append_build_context=true)
    output_format = documenter_output
    for NOTEBOOK_DIR in NOTEBOOK_DIRS
        bopts = BuildOptions(NOTEBOOK_DIR; output_format)
        build_notebooks(bopts, hopts)
    end
    return nothing
end

function mdfile(path...)
    joinpath(dirname(@__DIR__), "notebooks", path...)
end

# Build the notebooks; defaults to true.
if get(ENV, "BUILD_DOCS_NOTEBOOKS", "true") == "true"
    build()
end

sitename = "GiggleLiu's notebooks"
pages = [
    "PlutoStaticHTML" => "index.md",
    "Yao" => ["YaoBlocks"=>mdfile("yao, yaoblocks.md")],
    "Tensor Network" => [
            "Tropical Tensors"=>mdfile("tropical", "tropicaltensornetwork.md"),
            "Tropical GEMM"=>mdfile("tropical, tropicalgemm.md"),
        ]
]

# Using MathJax3 since Pluto uses that engine too.
mathengine = MathJax3()
prettyurls = get(ENV, "CI", nothing) == "true"
format = HTML(; mathengine, prettyurls)
modules = [PlutoStaticHTML]
strict = true
checkdocs = :none
makedocs(; sitename, pages, format, modules, strict, checkdocs)

deploydocs(;
    branch="docs-output",
    devbranch="main",
    repo="github.com/GiggleLiu/notebooks.git",
    push_preview=false
)
