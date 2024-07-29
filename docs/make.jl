# Add the path to the src directory of your package
push!(LOAD_PATH, "../src/")

using Documenter
using TaskBasedProduction

makedocs(
    sitename = "TaskBasedProduction.jl",
    modules  = [TaskBasedProduction],
    pages = [
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Examples" => "examples.md",
        "FAQ" => "faq.md",
        "Contributing" => "contributing.md"
    ]
)

deploydocs(
    repo = "github.com/haanwinckel/TaskBasedProduction",
    branch = "gh-pages",  # The branch to deploy to
    devbranch = "main",   # The development branch
    target = "docs"       # The directory to deploy to
)