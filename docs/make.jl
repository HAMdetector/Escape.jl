using Documenter
using Escape

makedocs(
    sitename = "Escape",
    format = Documenter.HTML(),
    modules = [Escape]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
