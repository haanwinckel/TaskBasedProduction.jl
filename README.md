# TaskBasedProduction

[![Build Status](https://github.com/lucalore98/TaskBasedProduction.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/lucalore98/TaskBasedProduction.jl/actions/workflows/CI.yml?query=branch%3Amain)

TaskBasedProduction is a Julia package that provides functions for calculating unit labor demands, marginal products of labor, production functions, assignment thresholds, elasticities of substitution, and complementarities among worker types in a task-based production model. The package includes utilities for handling incomplete gamma functions and power series representations to facilitate these calculations. This package was developed by Daniel Haanwinckel and Luca Lorenzini based on the paper *Supply, Demand, Institutions, and Firms: A Theory of Labor Market Sorting and the Wage Distribution* by Daniel Haanwinckel ([NBER Working Paper No. 31318](https://www.nber.org/papers/w31318)).

## Installation

To install TaskBasedProduction, you can clone the repository and add it to your Julia environment:

```julia
using Pkg
Pkg.add(url="https://github.com/haanwinckel/TaskBasedProduction.git")

using Revise
using TaskBasedProduction 
using SpecialFunctions
using QuadGK

θ = 1.0
κ = 0.5
z = 1.2
αVec = [0.1, 0.2, 0.3]
labor_input = [0.5; 0.04; 0.19]

initial_guess = find_initial_guess(θ, κ, z, αVec; threshold=1e-2)
q, xT, fval = prod_fun(labor_input, θ, κ, z, αVec; initial_guess=initial_guess, x_tol=1e-10)

println("Quantity Produced: ", q)
println("Task Thresholds: ", xT)
println("Approximation error: ", fval)

# Call unitInputDemand and print the output
labor_input2 = q * unitInputDemand(xT, θ, κ, z, αVec)
println("Labor Demand: ", labor_input2)
println("Error: ", labor_input2 - labor_input)

# Call margProdLabor with labor demand and print the output
MPL = margProdLabor(labor_input, θ, κ, z, αVec)
println("Marginal Products of Labor (with labor demand): ", MPL)

MPL = margProdLabor(labor_input, θ, κ, z, αVec, xT)
println("Marginal Products of Labor (with labor demand): ", MPL)

# Call elasticity_sub_comp with labor demand, MPL, xT, and parameters of the gamma function
ϵ_sub, ϵ_compl = elasticity_sub_comp(labor_input, θ, κ, z, αVec, MPL, xT)
println("Allen partial elasticity of substitution: ", ϵ_sub)
println("Hicks partial elasticity of substitution: ", ϵ_compl)

## GENERAL PARAMETERIZATION OF FUNCTIONS

# Define the density function b_g(x)
b_g(x) = (x^(κ-1) * exp(-x/θ)) / (θ^κ * gamma(κ))
e_h1(x) = exp(0.1*x)
e_h2(x) = exp(0.2*x)
e_h3(x) = exp(0.3*x)
e_h = [e_h1, e_h2, e_h3]  # Example e_h functions

initial_guess_gen = find_initial_guess_gen(z, b_g, e_h; threshold=1e-2, verbose=false)
q_gen, xT_gen, fval = prod_fun_general(labor_input, z, b_g, e_h; initial_guess=initial_guess_gen)

labor_input_general = q_gen * unitInputDemand_general(xT_gen, z, b_g, e_h)
println("Labor Demand: ", labor_input_general)
println("Is approximation close? ", isapprox(labor_input, labor_input_general, atol=1e-6))

println("Quantity Produced: ", q_gen)
println("Task Thresholds: ", xT_gen)
MPL_gen = margProdLabor_general(labor_input_general, z, b_g, e_h, xT_gen, q_gen)
ϵ_sub_gen, ϵ_compl_gen = elasticity_sub_comp_general(labor_input_general, z, b_g, e_h, MPL_gen, xT_gen)

println("General case Allen partial elasticity of substitution: ", ϵ_sub_gen)
println("General case Hicks partial elasticity of substitution: ", ϵ_compl_gen)
```
## Functions and Features
# 1) unitInputDemand
Calculates unit labor demands given blueprint scale θ, blueprint shape κ, productivity z, an array of comparative advantage values αVec, and an array xT of thresholds in task space.
```julia
 unitInputDemand(θ::Real, κ::Real, z::Real, αVec::Array{<:Real}, xT::Array{<:Real}, skipParamChecks::Bool = false)
 ```
# 2) margProdLabor
Calculates marginal products of labor for each worker type given blueprint characteristics or an array of labor demands.
```julia 
margProdLabor(inputDemand::Array{<:Real}, αVec::Array{<:Real}, xT::Array{<:Real})
```
# 3) prod_fun
Calculates the quantity produced (q) and task thresholds (xbar) given labor inputs (l), blueprint scale θ, blueprint shape κ, productivity z, and an array of comparative advantage values αVec with H elements (one for each worker type).

``` julia
prod_fun(l::AbstractArray{<:Real}, θ::Real, κ::Real, z::Real, αVec::AbstractArray{<:Real})
```
# 4) elasticity_sub_comp
Calculates the elasticity of substitution and complementarities for given labor inputs (l), an array xT of thresholds in task space (dimension H-1), array of marginal product of labor (MPL) for each H labor types, blueprint scale θ, blueprint shape κ, productivity z, and an array of comparative advantage values αVec with H elements (one for each worker type). The function returns two matrices representing the elasticity of substitution and complementarity values for each worker type h (rows) relative to worker type h_prime (columns).
``` julia
elasticity_sub_comp(xT::AbstractArray{<:Real}, l::AbstractArray{<:Real}, q::Real, MPL::AbstractArray{<:Real}, θ::Union{Real, AbstractArray{<:Real}}, κ::Union{Real, AbstractArray{<:Real}}, z::Union{Real, AbstractArray{<:Real}}, αVec::AbstractArray{<:Real})
```
# 5) unitInputDemand_general
Calculates unit labor demands given an array xT of H-1 thresholds in task space, productivity value z, a density function b_g for the task distribution, and an array e_h of H functions representing the cost of each labor type as a function of task complexity.

``` julia 
unitInputDemand_general(xT::Vector{Float64}, z::Real, b_g::Function, e_h::Vector{Function})
```
# 6) margProdLabor_general
Calculate marginal products of labor for each worker type given an array of H labor demand values, a vector of comparative advantage functions e_h, and an array of H-1 task thresholds xT that corresponds to that labor demand.
``` julia
margProdLabor_general(inputDemand::Array{<:Real}, αVec::Array{<:Real}, xT::Array{<:Real})
```
# 7) prod_fun_general
Calculates the quantity produced (q), and task thresholds (xbar) given labor inputs (l), productivity z, and general blueprint density function and a vector of efficiency functions, one for each labor type.

``` julia
prod_fun_general(l::AbstractArray{<:Real}, z::Real, b_g::Function, e_h::Vector{Function})
```
# 8) elasticity_sub_comp_general
Calculates the elasticity of substitution and complementarities for given labor inputs (l), an array xT of thresholds in task space (dimension H-1), array of marginal product of labor (MPL) for each H labor types, blueprint density function b_g, firm productivity z, and a vector of comparative advantage function e_h. The function returns two matrices representing the elasticity of substitution and complementarity values for each worker type h (rows) relative to worker type h_prime (columns).
``` julia
elasticity_sub_comp_general(xT::AbstractArray{<:Real}, l::AbstractArray{<:Real}, q::Real, MPL::AbstractArray{<:Real}, z::Real, b_g::Function, e_h::Vector{Function})
```

## Contributing
Contributions are welcome! Please feel free to submit a pull request or open an issue if you have any suggestions or find any bugs.

## License
This project is licensed under the MIT License - see the LICENSE file for details.