using Documenter
using Escape

makedocs(
    sitename = "Escape.jl",
    format = Documenter.HTML(
        canonical="https://HAMdetector.github.io/MyPkg2.jl",
    ),
    modules = [Escape],
    repo = "github.com/HAMdetector/Escape.jl.git",
    pages = [
        "quick_start.md",
        "advanced_usage.md"
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/HAMdetector/Escape.jl.git"
)
