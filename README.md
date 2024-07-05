# TaskBasedProduction

[![Build Status](https://github.com/lucalore98/TaskBasedProduction.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/lucalore98/TaskBasedProduction.jl/actions/workflows/CI.yml?query=branch%3Amain)

TaskBasedProduction is a Julia package that provides functions for calculating unit labor demands, marginal products of labor, production function, assignment thresholds, elasticities of substitution and complementarities among worker types in a task-based production model. The package includes utilities for handling incomplete gamma functions and power series representations to facilitate these calculations.

## Installation

To install TaskBasedProduction, you can clone the repository and add it to your Julia environment:

```julia
using Pkg
Pkg.add(url="https://github.com/lucalore98/TaskBasedProduction.jl.git")

## Usage Example

using Revise
using TaskBasedProduction 
using SpecialFunctions
using QuadGK
θ = 1.0
κ = 0.5
z = 1.2
αVec = [0.1, 0.2, 0.3]
xT = [0.4, 0.5]

# Call unitInputDemand and print the output
labor_demand = unitInputDemand(θ, κ, z, αVec, xT)
println("Labor Demand: ", labor_demand)

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
# Call elasticity_substitution with labor demand, MPL, xT and parameters of the gamma function
ϵ_sub, ϵ_compl=elasticity_sub_comp(xT, labor_demand, q, MPL, θ, κ, z, αVec)
println("Allen partial elasticity of substitution:", ϵ_sub)
println("Hicks partial elasticity of substitution:", ϵ_compl)

## GENERAL PARAMETERIZATION OF FUNCTIONS
# Example usage


# Define the density function b_g(x)
b_g(x) = (x^(κ-1) * exp(-x/θ)) / (θ^κ * gamma(κ))
e_h1(x)=exp(0.1*x)
e_h2(x)=exp(0.2*x)
e_h3(x)=exp(0.3*x)
e_h = [e_h1, e_h2, e_h3]  # Example e_h functions

labor_demand_general = unitInputDemand_general(xT, z,b_g, e_h)
println("Labor Demand: ", labor_demand_general)
isapprox(labor_demand, labor_demand_general, atol=1e-6)

q_gen, xT_gen= prod_fun_general(labor_demand_general,z,b_g, e_h)
println("Quantity Produced: ", q_gen)
println("Task Thresholds: ", xT_gen)
MPL_gen=margProdLabor_general(xT_gen, labor_demand_general, e_h)
ϵ_sub_gen, ϵ_compl_gen=elasticity_sub_comp_general(xT_gen, labor_demand_general, q_gen, MPL,z,b_g,e_h)
println("MPL General:", MPL_gen)
println("Allen partial elasticity of substitution, general case:", ϵ_sub_gen)
println("Hicks partial elasticity of substitution, general case:", ϵ_compl_gen)
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

5) unitInputDemand_general
Calculates unit labor demands given an array xT of H-1 thresholds in task space, productivity value z, 
a density function b_g for the task distribution, and an array e_h of H functions
representing the cost of each labor type as a function of task complexity.
unitInputDemand_general(xT::Vector{Float64}, z::Real, b_g:: Function, e_h::Vector{Function})

6)  margProdLabor_general
Calculate marginal products of labor for each worker type given an array
of H labor demand values, a vector of comparative advantage functions e_h, and
an array of H-1 task thresholds xT that corresponds to that labor demand.
margProdLabor_general(inputDemand::Array{<:Real}, αVec::Array{<:Real}, xT::Array{<:Real})

7) prod_fun_general

Calculates the quantity produced (q), and task thresholds (xbar)
given labor inputs (l),  productivity z, and general blueprint density function and a vector of efficiency functions, one for each labor type.
prod_fun_general(l::AbstractArray{<:Real}, z::Real, b_g:: Function, e_h::Vector{Function})

8) elasticity_sub_comp_general
Calculates the elasticity of substitution and complementarities for a given labor inputs (l), an array xT of thresholds in task space (dimension H-1), array of marginal product of labor (MPL) for each H labor types, blueprint density function b_g, firm productivity z, and a vector of comparative advantage function e_h. The function returns two matrices representing the elasticity of substitution and complementarity values for each worker type h (rows) relative to worker type h_prime (columns).
elasticity_sub_comp_general(xT::AbstractArray{<:Real},l::AbstractArray{<:Real}, q::Real, MPL::AbstractArray{<:Real}, z:: Real, b_g::Function, e_h::Vector{Function})

## Contributing
Contributions are welcome! Please feel free to submit a pull request or open an issue if you have any suggestions or find any bugs.

## License
This project is licensed under the MIT License - see the LICENSE file for details.


