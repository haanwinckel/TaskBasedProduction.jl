# Add the path to the src directory of your package
push!(LOAD_PATH, "../src/")

using Documenter
using TaskBasedProduction

makedocs(
    sitename = "TaskBasedProduction.jl",
    modules  = [TaskBasedProduction],
    pages = [
        "Home" => "index.md"
    ]
)
