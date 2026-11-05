# Documentation site of GarsideHomology.jl, built with Documenter.jl.
#
# Local build, from the root of the repository:
#     julia --project=docs -e "using Pkg; Pkg.instantiate()"
#     julia --project=docs docs/make.jl
# then open docs/build/index.html. Outside the CI, deploydocs publishes nothing; the site is
# published to GitHub Pages by the workflow .github/workflows/Documentation.yml.

using Documenter
using GarsideHomology
using Chevie   # loads the extension GarsideHomologyChevieExt, documented on its own page

# The tables of the pages Examples, generated from test/data/expected/*.toml.
include(joinpath(@__DIR__, "examples_tables.jl"))

# The extension module, under a name that the page reference_manual/chevie.md can use.
const GarsideHomologyChevieExt = Base.get_extension(GarsideHomology, :GarsideHomologyChevieExt)

makedocs(;
    sitename = "GarsideHomology.jl",
    authors = "Yves Guiraud",
    # The package with its submodule GarsideHomology.Complexes, and the extension.
    modules = [GarsideHomology, GarsideHomologyChevieExt],
    # Every docstring of an exported or public name must be in the manual; the docstrings of
    # the internal functions stay in the code. A missing docstring or a broken cross-reference
    # makes the build fail (warnonly = false, the default).
    checkdocs = :public,
    doctest = false,
    format = Documenter.HTML(;
        # Pretty URLs (germs/ instead of germs.html) on the CI only, so that the local build
        # can be read by opening docs/build/index.html, without a web server.
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://yvesguiraud.github.io/GarsideHomology.jl/stable/",
        edit_link = "main",
        # Hides the headings of the current page in the sidebar (docs/src/assets/custom.css).
        assets = ["assets/custom.css"],
    ),
    pages = [
        "Home" => "index.md",
        "User guide" => [
            "Getting started" => "user_guide/getting_started.md",
            "Background and conventions" => "user_guide/background.md",
            "Resolutions, cells and elements" => "user_guide/resolutions.md",
            "Display settings" => "user_guide/display.md",
            "Coefficients and homology" => "user_guide/homology.md",
            "Orderings, certificates and snapshots" => "user_guide/orderings.md",
            "Working with other packages" => "user_guide/interoperability.md",
        ],
        "Examples" => [
            "Overview" => "examples/overview.md",
            "Exceptional braid groups" => "examples/exceptional.md",
            "Spherical Artin groups" => "examples/spherical.md",
            "Imprimitive braid groups" => "examples/imprimitive.md",
        ],
        "Reference manual" => [
            "Overview" => "reference_manual/overview.md",
            "Germs" => "reference_manual/germs.md",
            "Orderings" => "reference_manual/orderings.md",
            "Cells" => "reference_manual/cells.md",
            "Resolutions and certificates" => "reference_manual/resolutions.md",
            "Coefficients and homology" => "reference_manual/coefficients.md",
            "Complexes" => "reference_manual/complexes.md",
            "Chevie extension" => "reference_manual/chevie.md",
        ],
        "Bibliography" => "bibliography.md",
    ],
)

# main gives dev/, a tag vX.Y.Z gives vX.Y.Z/ and stable/; nothing is published for a pull
# request. Authentication: the deploy key DOCUMENTER_KEY if the secret is set, otherwise
# GITHUB_TOKEN.
deploydocs(;
    repo = "github.com/yvesguiraud/GarsideHomology.jl.git",
    devbranch = "main",
)
