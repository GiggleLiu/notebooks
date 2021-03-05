using Documenter

# Pluto pages
import Pkg

Pkg.add([
Pkg.PackageSpec(url="https://github.com/GiggleLiu/PlutoUtils.jl", rev="static-export"),
Pkg.PackageSpec(url="https://github.com/fonsp/Pluto.jl", rev="05e5b68"),
]);

makedocs(;
    modules=Module[],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/GiggleLiu/NiLang.jl/blob/{commit}{path}#L{line}",
    sitename="Notebooks",
    authors="JinGuo Liu, thautwarm",
)

import PlutoUtils

PlutoUtils.Export.github_action(; notebook_dir=joinpath(dirname(@__DIR__), "notebooks"), offer_binder=false, export_dir=joinpath(@__DIR__, "build", "notebooks"), generate_default_index=false, project=@__DIR__)


deploydocs(;
    repo="github.com/GiggleLiu/notebooks.git",
)
