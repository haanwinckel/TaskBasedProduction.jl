using Revise
using TaskBasedProduction 
using SpecialFunctions
using QuadGK
θ = 1.0
κ = 0.5
z = 1.2
αVec = [0.1, 0.2, 0.3]
labor_input=[0.5178877620105745; 0.04097664559663147; 0.18579978864645658;;]

q, xT, fval= prod_fun(labor_input, θ, κ, z, αVec)
println("Quantity Produced: ", q)
println("Task Thresholds: ", xT)

# Call unitInputDemand and print the output
labor_input2 = unitInputDemand( xT, θ, κ, z, αVec)
println("Labor Demand: ", labor_input2)
labor_input2-labor_input

# Call margProdLabor with labor demand and print the output
MPL= margProdLabor(labor_input,  θ, κ, z, αVec)
println("Marginal Products of Labor (with labor demand): ",MPL)

# Call elasticity_substitution with labor demand, MPL, xT and parameters of the gamma function
ϵ_sub, ϵ_compl=elasticity_sub_comp(labor_input, θ, κ, z, αVec)
println("Allen partial elasticity of substitution:", ϵ_sub)
println("Hicks partial elasticity of substitution:", ϵ_compl)

## GENERAL PARAMETERIZATION OF FUNCTIONS

# Define the density function b_g(x)
b_g(x) = (x^(κ-1) * exp(-x/θ)) / (θ^κ * gamma(κ))
e_h1(x)=exp(0.1*x)
e_h2(x)=exp(0.2*x)
e_h3(x)=exp(0.3*x)
e_h = [e_h1, e_h2, e_h3]  # Example e_h functions
q_gen, xT_gen= prod_fun_general(labor_input,z,b_g, e_h)

labor_input_general = unitInputDemand_general(xT_gen, z, b_g, e_h)
println("Labor Demand: ", labor_input_general)
isapprox(labor_input, labor_input_general, atol=1e-6)

println("Quantity Produced: ", q_gen)
println("Task Thresholds: ", xT_gen)
MPL_gen=margProdLabor_general(labor_input_general, z, b_g, e_h)
ϵ_sub_gen, ϵ_compl_gen=elasticity_sub_comp_general(labor_input_general, z, b_g, e_h)



# This is how I would do this: unique input is labor_input 
MPL=margProdLabor(labor_input, θ, κ, z, αVec)

q, xT= prod_fun(labor_input, θ, κ, z, αVec)

qp, xTp= prod_fun(labor_input+1e-5*[1;0;0], θ, κ, z, αVec)

MPLp=margProdLabor(labor_input+1e-5*[1;0;0], θ, κ, z, αVec)

(qp-q)*1e5
(MPLp-MPL)*1e5

f13=ans[3]
f13*q/(MPL[1]*MPL[3])


function margProdLabor_old(inputDemand::AbstractArray{<:Real}, αVec::AbstractArray{<:Real}, xT::AbstractArray{<:Real})
    mpl_over_mpl1 = [1; cumprod(exp.(diff(αVec, dims=1) .* xT), dims=1)]
    mpl1 = 1 / sum(mpl_over_mpl1 .* inputDemand)
    mpl = mpl_over_mpl1 * mpl1
    return mpl, inputDemand
end

function margProdLabor_old(θ::Real, κ::Real, z::Real, αVec::Array{<:Real}, xT::Array{<:Real})
    l = unitInputDemand(xT, θ, κ, z, αVec)
    return margProdLabor_old(l, αVec, xT)
end

q2, xT2= prod_fun(labor_input, θ, κ, z, αVec)
MPL2, l2 = margProdLabor_old(θ, κ, z, αVec, xT2)

q-q2
xT-xT2
MPL-MPL2
labor_input-l2



#Now let's do the same after increasing employment of type 1 by 0.00001:
qp2, xTp2= prod_fun(labor_input+1e-5*[1;0;0], θ, κ, z, αVec)
MPLp2, lp2 = margProdLabor_old(θ, κ, z, αVec, xTp2)
#Below, I check that the change in output is approximately the same as would be predicted by MPL[1]:
(qp2-q2)*1e5
xTp-xTp2

#Now let's find, numerically, how much MPL[3] changes with the addition of type one workers: that is, the second derivative of the production function
(MPLp2-MPL2)*1e5
f13=ans[3]
#Let's use this second derivative to calculate the numerical elasticity of complementarity between types 1 and 3:
f13*q2/(MPL2[1]*MPL[3])

labor_perturbated=labor_input+1e-5*[1;0;0]
labor_perturbated-lp2


