using Revise
using TaskBasedProduction 
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