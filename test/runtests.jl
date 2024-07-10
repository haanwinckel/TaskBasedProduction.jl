using Test
using Revise
using TaskBasedProduction
using SpecialFunctions
# Hypothetical test for unitInputDemand
@testset "unitInputDemand Tests" begin
    θ = 1.0
    κ = 0.5
    z = 1.2
    αVec = [0.1, 0.2, 0.3]
    xT = [0.4, 0.5]

    expected_labor_demand = [0.5178877620105745; 0.04097664559663147; 0.18579978864645658;;]  # Replace with actual expected values
    labor_demand = unitInputDemand(θ, κ, z, αVec, xT)

    @test isapprox(labor_demand, expected_labor_demand, atol=1e-5)
end




# Test for margProdLabor with labor demand
@testset "margProdLabor with labor demand Tests" begin
    labor_demand =  [0.5178877620105745; 0.04097664559663147; 0.18579978864645658;;]
    αVec = [0.1, 0.2, 0.3]
    xT = [0.4, 0.5]

    expected_marginal_products = [1.309184899610229, 1.3626137489243066, 1.4324764497687001]  # Replace with actual expected values
    marginal_products = margProdLabor(labor_demand, αVec, xT)

    @test isapprox(marginal_products, expected_marginal_products, atol=1e-5)
end

# Test for margProdLabor with blueprint characteristics
@testset "margProdLabor with blueprint characteristics Tests" begin
    θ = 1.0
    κ = 0.5
    z = 1.2
    αVec = [0.1, 0.2, 0.3]
    xT = [0.4, 0.5]

    expected_marginal_products_alt = [1.309184899610229, 1.3626137489243066, 1.4324764497687001]  # Replace with actual expected values
    marginal_products_alt = margProdLabor(θ, κ, z, αVec, xT)

    @test isapprox(marginal_products_alt, expected_marginal_products_alt, atol=1e-5)
end




@testset "unitInputDemand Comparison Tests, Elasticity of complementarity" begin
    θ = 1.0
    κ = 0.5
    z = 1.2
    αVec = [0.1, 0.2, 0.3]
    xT = [0.4, 0.5]

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
    labor_demand_specific = unitInputDemand(θ, κ, z, αVec, xT)
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
    labor_demand_specific = unitInputDemand(θ, κ, z, αVec, xT)
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

