using Test
using Revise
using TaskBasedProduction
using SpecialFunctions
# Test for unitInputDemand
@testset "unitInputDemand Tests" begin
    θ = 1.0
    κ = 0.5
    z = 1.2
    αVec = [0.1, 0.2, 0.3]
    labor_input=[0.5178877620105745; 0.04097664559663147; 0.18579978864645658;;]
    q, xT= prod_fun(labor_input, θ, κ, z, αVec)  
    labor_input2 = unitInputDemand( xT, θ, κ, z, αVec)

    @test isapprox(labor_input, labor_input2, atol=1e-5)
end



@testset "unitInputDemand Comparison Tests, Elasticity of complementarity" begin
    θ = 1.0
    κ = 0.5
    z = 1.2
    αVec = [0.1, 0.2, 0.3]
    labor_input=[0.5178877620105745; 0.04097664559663147; 0.18579978864645658;;]

    # Define the density function b_g(x)
    b_g(x) = (x^(κ-1) * exp(-x/θ)) / (θ^κ * gamma(κ))
   

    # Define the e_h functions
    e_h1(x) =exp(αVec[1]* x)
    e_h2(x) = exp(αVec[2] * x)
    e_h3(x) = exp(αVec[3] * x)
    e_h = [e_h1, e_h2, e_h3]
    q_gen, xT_gen= prod_fun_general(labor_input,z,b_g, e_h)
    q, xT= prod_fun(labor_input, θ, κ, z, αVec)
    # Test to check whether q is the same 
    @test isapprox(q, q_gen, atol=1e-5)
    # Test to check whether xT is the same 
    @test isapprox(xT, xT_gen, atol=1e-5)
    # Compute labor demand using the general function
    labor_input_general = unitInputDemand_general(xT_gen, z, b_g, e_h)
    # Compute labor demand using the specific function
    labor_input_specific = unitInputDemand( xT, θ, κ, z, αVec)
    # Test to check if the two labor inputs are the same
    @test isapprox(labor_input_general, labor_input_specific, atol=1e-5)
    
    MPL= margProdLabor(labor_input_specific,  θ, κ, z, αVec)
    MPL_gen=margProdLabor_general(labor_input_general, z, b_g, e_h)
    # Test the two MPLs to be the same
    @test isapprox(MPL,MPL_gen,  atol=1e-4)
    ϵ_sub, ϵ_compl=elasticity_sub_comp(labor_input_specific, θ, κ, z, αVec)
    ϵ_sub_gen, ϵ_compl_gen=elasticity_sub_comp_general(labor_input_general, z, b_g, e_h)
    @test isapprox(ϵ_sub,ϵ_sub_gen, atol=1e-4)
    @test isapprox(ϵ_compl,ϵ_compl_gen, atol=1e-4)
    
   
    # Test to check if an error is thrown when b_g is not a valid density function
    b_g_invalid = x -> 4*exp(-x)  # This does not integrate to 1 over the entire domain
    @test_throws ErrorException unitInputDemand_general(xT_gen, z, b_g, e_h) 
    
    function numerical_second_deriv(labor_input, θ::Real, κ::Real, z::Real, αVec::AbstractVector{<:Real}, h::Int, hprime::Int; hstep=1e-5)
        @assert h > 0 "h must be a natural number (positive integer)"
        @assert hprime > 0 "hprime must be a natural number (positive integer)"
        
        perturbation = zeros(length(labor_input))
        perturbation[h] = hstep
        
        MPL_plus_h = margProdLabor(labor_input + perturbation, θ, κ, z, αVec)
        MPL_minus_h = margProdLabor(labor_input - perturbation, θ, κ, z, αVec)
        
        return (MPL_plus_h[hprime] - MPL_minus_h[hprime]) / (2 * hstep)
    end
    
    function numerical_elasticity_compl(labor_input, θ::Real, κ::Real, z::Real, αVec::AbstractVector{<:Real})
        MPL = margProdLabor(labor_input, θ, κ, z, αVec)
        H = length(MPL)
        ϵ_compl_numerical = zeros(Float64, H, H)
        
        q, xT= prod_fun(labor_input, θ, κ, z, αVec)  
        
        for h in 1:H
            for hprime in h+1:H
                MPL_d = numerical_second_deriv(labor_input, θ, κ, z, αVec, h, hprime)
                ϵ_compl_numerical[h, hprime] = q * MPL_d / (MPL[h] * MPL[hprime])
            end
        end
        
        return ϵ_compl_numerical
    end
    

    ϵ_compl2 = numerical_elasticity_compl(labor_input_specific, θ, κ, z, αVec)
    println(ϵ_compl2)
    ϵ_compl2=numerical_elasticity_compl(labor_input_specific, θ, κ, z, αVec)
    @test isapprox(ϵ_compl,ϵ_compl2, atol=1e-3)
    @test isapprox(ϵ_compl_gen,ϵ_compl2, atol=1e-3)
end
#Different labor demand
@testset "unitInputDemand Comparison Tests, Elasticity of complementarity" begin
    θ = 1.0
    κ = 0.5
    z = 1.2
    αVec = [0.1, 0.2, 0.3]
    xT = [1.4, 1.5]

    # Define the density function b_g(x)
    b_g(x) = (x^(κ-1) * exp(-x/θ)) / (θ^κ * gamma(κ))
   

    # Define the e_h functions
    e_h1(x) =exp(αVec[1]* x)
    e_h2(x) = exp(αVec[2] * x)
    e_h3(x) = exp(αVec[3] * x)
    e_h = [e_h1, e_h2, e_h3]

    # Compute labor demand using the general function
    labor_demand_general = unitInputDemand_general(xT, z, b_g, e_h)
    # Compute labor demand using the specific function
    labor_demand_specific = unitInputDemand( xT, θ, κ, z, αVec)
    MPL= margProdLabor(labor_demand_specific, αVec, xT)
    q, xT= prod_fun(labor_demand_specific, θ, κ, z, αVec)
    q_gen, xT_gen= prod_fun_general(labor_demand_general,z,b_g, e_h)
    MPL_gen=margProdLabor_general(xT_gen, labor_demand_general, e_h)
    ϵ_sub, ϵ_compl=elasticity_sub_comp(xT, labor_demand_specific, q, MPL, θ, κ, z, αVec)
    ϵ_sub_gen, ϵ_compl_gen=elasticity_sub_comp_general(xT_gen, labor_demand_general, q_gen, MPL_gen,z,b_g,e_h)
    @test isapprox(MPL,MPL_gen,  atol=1e-4)
    @test isapprox(ϵ_sub,ϵ_sub_gen, atol=1e-4)
    @test isapprox(ϵ_compl,ϵ_compl_gen, atol=1e-4)
    # Test to check if the two outputs are the same
    @test isapprox(labor_demand_general, labor_demand_specific, atol=1e-5)
    # Test to check if an error is thrown when b_g is not a valid density function
    b_g_invalid = x -> 4*exp(-x)  # This does not integrate to 1 over the entire domain
    @test_throws ErrorException unitInputDemand_general(xT, z, b_g_invalid, e_h)
   
    
    # Test to check whether q is the same 
    @test isapprox(q, q_gen, atol=1e-5)
    # Test to check whether xT is the same 
    @test isapprox(xT, xT_gen, atol=1e-5)


    function numerical_second_deriv(labor_demand_specific, margProdLabor, θ, κ, z, αVec, xT, h, hprime; hstep=1e-5)
        perturbation = zeros(length(xT)+1)
        perturbation[h] = hstep
        qp, xTp= prod_fun(labor_demand_specific+perturbation, θ, κ, z, αVec)
        qpp, xTpp= prod_fun(labor_demand_specific-perturbation, θ, κ, z, αVec)
              
        MPL_plus_h = margProdLabor(θ, κ, z, αVec, xTp)
        MPL_minus_h = margProdLabor(θ, κ, z, αVec, xTpp)
        
        return (MPL_plus_h[hprime] - MPL_minus_h[hprime]) / (2 * hstep)
    end
    
    function numerical_elasticity_compl(q, MPL, numerical_second_deriv, θ, κ, z, αVec, xT)
        H = length(MPL)
        ϵ_compl_numerical = zeros(Float64, H, H)
        
        for h in 1:H
            for hprime in h+1:H
                if h < hprime
                    MPL_d = numerical_second_deriv(labor_demand_specific, margProdLabor, θ, κ, z, αVec, xT, h, hprime)
                    ϵ_compl_numerical[h, hprime] = q * MPL_d / (MPL[h] * MPL[hprime])
                else
                    ϵ_compl_numerical[h, hprime] = NaN
                end
            end
        end
        
        return ϵ_compl_numerical
    end
    ϵ_compl2=numerical_elasticity_compl(q, MPL, numerical_second_deriv, θ, κ, z, αVec, xT)
    @test isapprox(ϵ_compl,ϵ_compl2, atol=1e-3)
    @test isapprox(ϵ_compl_gen,ϵ_compl2, atol=1e-3)
end

