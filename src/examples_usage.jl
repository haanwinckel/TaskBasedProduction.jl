using TaskBasedProduction 
θ = 1.0
κ = 0.5
z = 1.2
αVec = [0.1, 0.2, 0.3]
labor_input=[0.5; 0.04; 0.19;;]

initial_guess=find_initial_guess(θ, κ, z, αVec; threshold=1e-2)
q, xT, fval = prod_fun(labor_input, θ, κ, z, αVec; initial_guess=initial_guess,  x_tol=1e-10)

println("Quantity Produced: ", q)
println("Task Thresholds: ", xT)
println("Approximation error: ", fval)
# Call unitInputDemand and print the output
labor_input2 = unitInputDemand( xT, q, θ, κ, z, αVec)
println("Labor Demand: ", labor_input2)
println("Error", labor_input2-labor_input)

# Call margProdLabor with labor demand and print the output
MPL= margProdLabor(labor_input,  θ, κ, z, αVec)
println("Marginal Products of Labor (with labor demand): ",MPL)

MPL= margProdLabor(labor_input,  θ, κ, z, αVec, xT)
println("Marginal Products of Labor (with labor demand): ",MPL)
MPL= margProdLabor(nothing, θ, κ, z, αVec)
println("Marginal Products of Labor (with labor demand): ",MPL)

# Call elasticity_substitution with labor demand, MPL, xT and parameters of the gamma function
ϵ_sub, ϵ_compl=elasticity_sub_comp(labor_input, θ, κ, z, αVec, MPL, xT)
println("Allen partial elasticity of substitution:", ϵ_sub)
println("Hicks partial elasticity of substitution:", ϵ_compl)

## GENERAL PARAMETERIZATION OF FUNCTIONS

# Define the density function b_g(x)
using SpecialFunctions # Note: this package is needed only to compare the output of the general parameterization with the output of the Gamma parameterization. It is used in the next line to obtain the PDF of a Gamma distribution.

b_g(x) = (x^(κ - 1) * exp(-x / θ)) / (θ^κ * gamma(κ)) # Gamma PDF with the same parameterization as above for comparability
e_h1(x)=exp(0.1*x)
e_h2(x)=exp(0.2*x)
e_h3(x)=exp(0.3*x)
e_h = [e_h1, e_h2, e_h3]  # Example e_h functions

initial_guess_gen=find_initial_guess_gen(z, b_g, e_h; threshold=1e-2, verbose=false)
q_gen, xT_gen,fval= prod_fun_general(labor_input,z,b_g, e_h; initial_guess=initial_guess_gen)

labor_input_general = unitInputDemand_general(xT_gen, q_gen, z, b_g, e_h)
println("Labor Demand: ", labor_input_general)
isapprox(labor_input, labor_input_general, atol=1e-6)

println("Quantity Produced: ", q_gen)
println("Task Thresholds: ", xT_gen)
MPL_gen=margProdLabor_general(labor_input_general, z, b_g, e_h, xT_gen, q_gen)
ϵ_sub_gen, ϵ_compl_gen=elasticity_sub_comp_general(labor_input_general, z, b_g, e_h, MPL_gen, xT_gen)



