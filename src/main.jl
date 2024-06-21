using TaskBasedProduction 
θ = 1.0
κ = 0.5
z = 1.2
αVec = [0.1, 0.2, 0.3]
xT = [0.4, 0.5]

# Call unitInputDemand and print the output
labor_demand = unitInputDemand(θ, κ, z, αVec, xT)
println("Labor Demand: ", labor_demand)

# Call component_positive_ups and print the output
Υ = 0.2
xT_low = 0.0
xT_high = 0.5
positive_ups = component_positive_ups(Υ, κ, xT_low, xT_high)
println("Component Positive UPS: ", positive_ups)

# Call component_negative_ups and print the output
negative_ups = component_negative_ups(Υ, κ, xT_low, xT_high)
println("Component Negative UPS: ", negative_ups)

# Call margProdLabor with labor demand and print the output
marginal_products = margProdLabor(labor_demand, αVec, xT)
println("Marginal Products of Labor (with labor demand): ", marginal_products)

# Call margProdLabor with blueprint characteristics and print the output
marginal_products_alt = margProdLabor(θ, κ, z, αVec, xT)
println("Marginal Products of Labor (with blueprint characteristics): ", marginal_products_alt)