using TaskBasedProduction 

θ = 1.0
κ = 0.5
z = 1.2
αVec = [0.1, 0.2, 0.3]
labor_input=[0.5; 0.04; 0.19;;]

initial_guess=find_initial_guess(θ, κ, z, αVec; threshold=1e-2)
q, xT = prod_fun(labor_input, θ, κ, z, αVec; initial_guess=initial_guess,  x_tol=1e-10)

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


# Call elasticity_substitution with labor demand, MPL, xT and parameters of the gamma function
ϵ_sub, ϵ_compl=elasticity_sub_comp(labor_input, θ, κ, z, αVec, MPL, xT)
println("Allen partial elasticity of substitution:", ϵ_sub)
println("Hicks partial elasticity of substitution:", ϵ_compl)
# Call elasticity_substitution with labor demand, MPL, xT and parameters of the gamma function
ϵ_sub, ϵ_compl=elasticity_sub_comp(labor_input, θ, κ, z, αVec, MPL, xT, q)
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
q_gen, xT_gen= prod_fun_general(labor_input,z,b_g, e_h; initial_guess=initial_guess_gen)

labor_input_general = unitInputDemand_general(xT_gen, q_gen, z, b_g, e_h)
println("Labor Demand: ", labor_input_general)
isapprox(labor_input, labor_input_general, atol=1e-6)

println("Quantity Produced: ", q_gen)
println("Task Thresholds: ", xT_gen)
MPL_gen=margProdLabor_general(labor_input_general, z, b_g, e_h, xT_gen, q_gen)
ϵ_sub_gen, ϵ_compl_gen=elasticity_sub_comp_general(nothing, z, b_g, e_h, MPL_gen, xT_gen, q_gen)

## Use case 1: competitive market with functional forms
q = 1
wage = [0.1; 0.2; 0.7]

# Compute thresholds xT suing the functional form 
diff_alpha = diff(αVec)
log_wage_ratios = log.(wage[2:end] ./ wage[1:end-1])
xT = (1 ./ diff_alpha) .* log_wage_ratios

# Calculate labor unit input requirements 
labor_input_1 = unitInputDemand(xT, q, θ, κ, z, αVec)
println("Labor Input: ", labor_input_1)
# Competitive labor market with general functions
# Objective function for optimization
function objective(x, h)
    x=x[1]
    return (e_h[h+1](x) / e_h[h](x) - wage[h+1] / wage[h])^2
end

# Find the solutions for xT_h
H=length(wage)
xT = Vector{Float64}(undef, H-1)
# Initial guess: recall that the first element is q
initial_guess=find_initial_guess_gen(z, b_g, e_h)
for h in 1:H-1
     x0=[initial_guess[h+1]]
    # Perform the optimization
    res = optimize(x ->objective(x, h), x0, LevenbergMarquardt())
    # Extract the optimized value
    xT[h]=res.minimizer[1]
end

labor_input_2 = unitInputDemand_general(xT, q, z, b_g, e_h)
println("Labor Input: ", labor_input_2)
# Use case 2: Given labor input, use the production function to obtain total production and task thresholds

labor_input=[0.5; 0.04; 0.19;;]
initial_guess=find_initial_guess(θ, κ, z, αVec; threshold=1e-2)
q, xT = prod_fun(labor_input, θ, κ, z, αVec; initial_guess=initial_guess)


#  General parameterization
initial_guess=find_initial_guess_gen(z, b_g, e_h)
q, xT=prod_fun_general(labor_input,z,b_g, e_h; initial_guess=initial_guess)

# Use case 3: Calculate the elasticity of complementarity and substitution given labor inputs and total production and task thresholds
# Call elasticity_substitution with labor demand given. 
ϵ_sub, ϵ_compl=elasticity_sub_comp(labor_input, θ, κ, z, αVec)
println("Allen partial elasticity of substitution:", ϵ_sub)
println("Hicks partial elasticity of substitution:", ϵ_compl)


# Use case 4: problem of the firm in a monopsonistic labor market 
# Define parameter
β=4
# Define general equilibrium objects
L=[1 ; 1 ; 1] # Total labor force 
p=1  # Price for the good
w_inclusive=[0.4; 0.9; 2]  # Inclusive value of wages 

# Define the objective function for optimization derived from the firm problem 
function objective_to_minimize(initial_guess)
    # Compute q and xT from the initial guess
    q = exp(initial_guess[1])
    xT = cumsum(exp.(initial_guess[2:end]))
    
    # Calculate labor input demand
    labor_input = unitInputDemand(xT, q, θ, κ, z, αVec)
    
    # Calculate MPL
    MPL = margProdLabor(labor_input, θ, κ, z, αVec, xT, q)
    
    # Calculate wages
    w = p * (β / (β + 1)) * MPL
    
    # Calculate labor supply
    labor_supply = (w ./ w_inclusive) .^ β .* L
    
    # Calculate error
    err = log.(labor_input ./ labor_supply)
    
    # Objective to minimize: sum of squared errors
    return sum(err .^ 2)
end

# Initial guess
initial_guess = find_initial_guess(θ, κ, z, αVec; threshold=1e-2)
using LeastSquaresOptim
result = optimize(objective_to_minimize, initial_guess, LevenbergMarquardt())

# Extract results
optimal_initial_guess = result.minimizer
q_opt = exp(optimal_initial_guess[1])
xT_opt = cumsum(exp.(optimal_initial_guess[2:end]))

# Print results
println("Optimal q: ", q_opt)
println("Optimal xT: ", xT_opt)

# Same thing but with the functions handling general blueprint functions and efficiency functions.
function objective_to_minimize(initial_guess)
    # Compute q and xT from the initial guess
    q = exp(initial_guess[1])
    xT = cumsum(exp.(initial_guess[2:end]))
    
    # Calculate labor input demand
    labor_input = unitInputDemand_general(xT, q, z, b_g, e_h)
    
    # Calculate MPL
    MPL= margProdLabor_general(labor_input, z, b_g, e_h, xT, q)
    
    # Calculate wages
    w = p * (β / (β + 1)) * MPL
    
    # Calculate labor supply
    labor_supply = (w ./ w_inclusive) .^ β .* L
    
    # Calculate error
    err = log.(labor_input ./ labor_supply)
    
    # Objective to minimize: sum of squared errors
    return sum(err .^ 2)
end

using LeastSquaresOptim
result = optimize(objective_to_minimize, initial_guess, LevenbergMarquardt())

# Extract results
optimal_initial_guess = result.minimizer
q_opt = exp(optimal_initial_guess[1])
xT_opt = cumsum(exp.(optimal_initial_guess[2:end]))

# Print results
println("Optimal q: ", q_opt)
println("Optimal xT: ", xT_opt)
