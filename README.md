# TaskBasedProduction

[![Build Status](https://github.com/lucalore98/TaskBasedProduction.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/lucalore98/TaskBasedProduction.jl/actions/workflows/CI.yml?query=branch%3Amain)

TaskBasedProduction is a Julia package that provides functions for calculating unit labor demands and marginal products of labor in a task-based production model. The package includes utilities for handling incomplete gamma functions and power series representations to facilitate these calculations.

## Installation

To install TaskBasedProduction, you can clone the repository and add it to your Julia environment:

```julia
using Pkg
Pkg.add(url="https://github.com/lucalore98/TaskBasedProduction.jl.git")

## Usage Example

using TaskBasedProduction

θ = 1.0
κ = 0.5
z = 1.2
αVec = [0.1, 0.2, 0.3]
xT = [0.4, 0.5]

# Call unitInputDemand and print the output
labor_demand = unitInputDemand(θ, κ, z, αVec, xT)
println("Labor Demand: ", labor_demand)
# Call prod_fun and print the output
q, xT= prod_fun(labor_demand, θ, κ, z, αVec)
println("Quantity Produced: ", q)
println("Task Thresholds: ", xT)
# Call margProdLabor with labor demand and print the output
MPL= margProdLabor(labor_demand, αVec, xT)
println("Marginal Products of Labor (with labor demand): ",MPL)

# Call margProdLabor with blueprint characteristics and print the output
MPL_alt = margProdLabor(θ, κ, z, αVec, xT)
println("Marginal Products of Labor (with blueprint characteristics): ", MPL_alt)
q=1
# Call elasticity_sub_comp with labor demand, MPL, xT and parameters of the gamma function and print the two outputs
ϵ_sub, ϵ_compl=elasticity_sub_comp(xT, labor_demand, q, MPL, θ, κ, z, αVec)
println("Allen partial elasticity of substitution:", ϵ_sub)
println("Hicks partial elasticity of substitution:", ϵ_compl)
```
## Functions and Features
1) unitInputDemand
Calculates unit labor demands given blueprint scale θ, blueprint shape κ, productivity z, an array of comparative advantage values αVec, and an array xT of thresholds in task space.
unitInputDemand(θ::Real, κ::Real, z::Real, αVec::Array{<:Real}, xT::Array{<:Real}, skipParamChecks::Bool = false)

2) margProdLabor
Calculates marginal products of labor for each worker type given blueprint characteristics or an array of labor demands. 
margProdLabor(inputDemand::Array{<:Real}, αVec::Array{<:Real}, xT::Array{<:Real})
margProdLabor(θ::Real, κ::Real, z::Real, αVec::Array{<:Real}, xT::Array{<:Real})

3) prod_fun
Calculates the quantity produced (q) and task thresholds (xbar) given labor inputs (l), blueprint scale θ, blueprint shape κ, productivity z, and an array of comparative advantage values αVec with H elements (one for each worker type).
 prod_fun(l::AbstractArray{<:Real}, θ::Real, κ::Real, z::Real, αVec::AbstractArray{<:Real})

4) elasticity_sub_comp
Calculates the elasticity of substitution and complementarities for a given labor inputs (l), an array xT of thresholds in task space (dimension H-1), array of marginal product of labor (MPL) for each H labor types, blueprint scale θ, blueprint shape κ, productivity z, and an array of comparative advantage values αVec with H elements (one for each worker type). Note that θ, κ, and z are not restricted to be scalar. The function returns two matrices representing the elasticity of substitution and complementarity values for each worker type h (rows) relative to worker type h_prime (columns).
elasticity_sub_comp(xT::AbstractArray{<:Real}, l::AbstractArray{<:Real}, q::Real, MPL::AbstractArray{<:Real}, θ::Union{Real, AbstractArray{<:Real}}, κ::Union{Real, AbstractArray{<:Real}}, z::Union{Real, AbstractArray{<:Real}}, αVec::AbstractArray{<:Real})


## Contributing
Contributions are welcome! Please feel free to submit a pull request or open an issue if you have any suggestions or find any bugs.

## License
This project is licensed under the MIT License - see the LICENSE file for details.


