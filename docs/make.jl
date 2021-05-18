using Documenter
using Escape

makedocs(
    sitename = "Escape",
    format = Documenter.HTML(),
    modules = [Escape],
    pages = [
        "installation.md",
        "quick_start.md",
        "advanced_usage.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
