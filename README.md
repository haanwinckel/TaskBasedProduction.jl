# TaskBasedProduction

[![Build Status](https://github.com/lucalore98/TaskBasedProduction.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/lucalore98/TaskBasedProduction.jl/actions/workflows/CI.yml?query=branch%3Amain)

TaskBasedProduction is a Julia package that provides functions for calculating unit labor demands and marginal products of labor in a task-based production model. The package includes utilities for handling incomplete gamma functions and power series representations to facilitate these calculations.

## Installation

To install TaskBasedProduction, you can clone the repository and add it to your Julia environment:

```julia
using Pkg
Pkg.add(url="https://github.com/lucalore98/TaskBasedProduction.jl.git")

## Usage
#Example: Calculating Unit Labor Demand

using TaskBasedProduction

# Define input parameters
θ = 1.0
κ = 0.5
z = 1.2
αVec = [0.1, 0.2, 0.3]
xT = [0.4, 0.5]

# Call unitInputDemand and print the output
labor_demand = unitInputDemand(θ, κ, z, αVec, xT)
println("Labor Demand: ", labor_demand)

#Example: Calculating Marginal Products of Labor
using TaskBasedProduction

# Define input parameters
θ = 1.0
κ = 0.5
z = 1.2
αVec = [0.1, 0.2, 0.3]
xT = [0.4, 0.5]

# Call margProdLabor with blueprint characteristics and print the output
marginal_products = margProdLabor(θ, κ, z, αVec, xT)
println("Marginal Products of Labor: ", marginal_products)
```
## Functions and Features
1) unitInputDemand
Calculates unit labor demands given blueprint scale θ, blueprint shape κ, productivity z, an array of comparative advantage values αVec, and an array xT of thresholds in task space.
unitInputDemand(θ::Real, κ::Real, z::Real, αVec::Array{<:Real}, xT::Array{<:Real}, skipParamChecks::Bool = false)

2) margProdLabor
Calculates marginal products of labor for each worker type given blueprint characteristics or an array of labor demands. 
margProdLabor(inputDemand::Array{<:Real}, αVec::Array{<:Real}, xT::Array{<:Real})
margProdLabor(θ::Real, κ::Real, z::Real, αVec::Array{<:Real}, xT::Array{<:Real})

3) component_positive_ups
Calculates the positive component of the UPS using the incomplete gamma function. This function is used internally by unitInputDemand.

4) component_negative_ups
Calculates the negative component of the UPS using the power series representation. This function is used internally by unitInputDemand.


## Contributing
Contributions are welcome! Please feel free to submit a pull request or open an issue if you have any suggestions or find any bugs.

## License
This project is licensed under the MIT License - see the LICENSE file for details.


