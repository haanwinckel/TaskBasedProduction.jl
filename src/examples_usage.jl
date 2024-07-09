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
isapprox(ϵ_sub,ϵ_sub_gen, atol=1e-4)
isapprox(ϵ_compl,ϵ_compl_gen, atol=1e-4)
isapprox(MPL,MPL_gen,  atol=1e-4)


# Function to compute the second derivative of the marginal products of labor
function compute_second_derivatives(labor_demand, MPL, h, hprime)
    delta=1e-12
    perturbed_quantities = copy(labor_demand)
    perturbed_quantities[hprime] += delta
    marginal_product_original =MPL[h]
    q2, xT_perturbated= prod_fun(perturbed_quantities, θ, κ, z, αVec)
    MPL_alt =margProdLabor(perturbed_quantities, αVec, xT_perturbated)
    marginal_product_perturbed=MPL_alt[h]
    return (marginal_product_perturbed - marginal_product_original) / delta
end

# Function to compute the elasticities of complementarity
function compute_elasticities(labor_demand, q, MPL)
    H = length(labor_demand)
    elasticities = zeros(H, H)
    for h in 1:H
        for hprime in 1:H
                second_derivative = compute_second_derivatives(labor_demand, MPL, h, hprime)
                if h<hprime 
                    elasticities[h, hprime] = q* second_derivative / (MPL[h] * MPL[hprime])
                else elasticities[h, hprime]=NaN
                end
            
        end
    end
    return elasticities
end

compute_elasticities(labor_demand, q, MPL)

h=2
hprime=3
delta=1e-5
perturbed_quantities = copy(labor_demand)
perturbed_quantities[hprime] += delta
q2, xT_perturbated= prod_fun(perturbed_quantities, θ, κ, z, αVec)
(q2-q)/delta
MPL_alt =margProdLabor(perturbed_quantities, αVec, xT_perturbated)
MPL_d=(MPL_alt[h]-MPL[h])/delta
q*MPL_d/(MPL[h]*MPL[hprime])