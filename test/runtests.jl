using Test
using TaskBasedProduction
using SpecialFunctions
# Function to run the tests
function run_tests(labor_input::Vector{Float64})
    # Initialize common parameters for tests
    θ = 1.0
    κ = 0.5
    z = 1.2
    αVec = [0.1, 0.2, 0.3]

    initial_guess = find_initial_guess(θ, κ, z, αVec; threshold=1e-2)
    q, xT, fval = prod_fun(labor_input, θ, κ, z, αVec; initial_guess=initial_guess, x_tol=1e-10)
    MPL = margProdLabor(labor_input, θ, κ, z, αVec)
    ϵ_sub, ϵ_compl = elasticity_sub_comp(labor_input, θ, κ, z, αVec, MPL, xT)

    # Test for unitInputDemand
    @testset "unitInputDemand Tests" begin
        labor_input2 = q * unitInputDemand(xT, θ, κ, z, αVec)
        @test isapprox(labor_input, labor_input2, atol=1e-5)
    end

    # Test for MPL function
    @testset "MPL function" begin
        MPL2 = margProdLabor(labor_input, θ, κ, z, αVec, xT)
        @test isapprox(MPL, MPL2, atol=1e-5)
    end

    # Numerical comparison for MPL
    @testset "MPL numerical comparison" begin
        H = length(labor_input)
        tol = 1e-5
        num_MPL = zeros(H)
        for i in 1:H
            perturbation = zeros(H)
            perturbation[i] = tol
            qp, xTp, fvalp = prod_fun(labor_input + perturbation, θ, κ, z, αVec; initial_guess=initial_guess, x_tol=1e-10)
            qn, xTn, fvaln = prod_fun(labor_input - perturbation, θ, κ, z, αVec; initial_guess=initial_guess, x_tol=1e-10)
            num_MPL[i] = (qp - qn) / (2 * tol)
        end
        @test isapprox(MPL, num_MPL, atol=1e-3)
    end

    # Define the density function b_g(x)
    b_g(x) = (x^(κ-1) * exp(-x/θ)) / (θ^κ * gamma(κ))
    e_h1(x) = exp(0.1 * x)
    e_h2(x) = exp(0.2 * x)
    e_h3(x) = exp(0.3 * x)
    e_h = [e_h1, e_h2, e_h3]  # Example e_h functions
    initial_guess_gen = find_initial_guess_gen(z, b_g, e_h; threshold=1e-2, verbose=false)
    q_gen, xT_gen, fval = prod_fun_general(labor_input, z, b_g, e_h; initial_guess=initial_guess_gen)
    labor_input_general = q_gen * unitInputDemand_general(xT_gen, z, b_g, e_h)
    MPL_gen = margProdLabor_general(labor_input_general, z, b_g, e_h, xT_gen, q_gen)
    ϵ_sub_gen, ϵ_compl_gen = elasticity_sub_comp_general(labor_input_general, z, b_g, e_h, MPL_gen, xT_gen)

    @testset "prod_fun_general comparison test" begin
        @test isapprox(q_gen, q, atol=1e-5)
        @test isapprox(xT_gen, xT, atol=1e-5)
    end

    @testset "Labor input general" begin
        @test isapprox(labor_input_general, labor_input, atol=1e-5)
    end

    @testset "MPL general" begin
        @test isapprox(MPL_gen, MPL, atol=1e-5)
    end

    @testset "Elasticity compl general" begin
        @test isapprox(ϵ_compl_gen, ϵ_compl, atol=1e-2)
    end

    @testset "Elasticity sub general" begin
        @test isapprox(ϵ_sub_gen, ϵ_sub, atol=10)
    end

    @testset "Invalid density" begin
        # Test to check if an error is thrown when b_g is not a valid density function
        b_g_invalid = x -> 4 * exp(-x)  # This does not integrate to 1 over the entire domain
        @test_throws ErrorException unitInputDemand_general(xT_gen, z, b_g_invalid, e_h)
    end

    @testset "Numerical Elasticity of complementarity" begin
        function numerical_second_deriv(labor_input, θ::Real, κ::Real, z::Real, αVec::AbstractVector{<:Real}, h::Int, hprime::Int; hstep=1e-5)
            @assert h > 0 "h must be a natural number (positive integer)"
            @assert hprime > 0 "hprime must be a natural number (positive integer)"

            perturbation = zeros(length(labor_input))
            perturbation[h] = hstep

            MPL_plus_h = margProdLabor(labor_input + perturbation, θ, κ, z, αVec)
            MPL_minus_h = margProdLabor(labor_input - perturbation, θ, κ, z, αVec)

            return (MPL_plus_h[hprime] - MPL_minus_h[hprime]) / (2 * hstep)
        end

        function numerical_elasticity_compl(labor_input, MPL, q, θ::Real, κ::Real, z::Real, αVec::AbstractVector{<:Real})
            H = length(MPL)
            ϵ_compl_numerical = zeros(Float64, H, H)
            for h in 1:H
                for hprime in h+1:H
                    MPL_d = numerical_second_deriv(labor_input, θ, κ, z, αVec, h, hprime)
                    ϵ_compl_numerical[h, hprime] = q * MPL_d / (MPL[h] * MPL[hprime])
                end
            end

            return ϵ_compl_numerical
        end

        ϵ_compl2 = numerical_elasticity_compl(labor_input, MPL, q, θ, κ, z, αVec)
        @test isapprox(ϵ_compl, ϵ_compl2, atol=1e-5)
        @test isapprox(ϵ_compl_gen, ϵ_compl2, atol=1e-5)
    end
end

# Initial labor inputs to test
labor_inputs = [[0.5; 0.04; 0.19], [0.2; 0.05; 0.5]]

# Run tests for each initial labor input
for labor_input in labor_inputs
    println("Running tests for labor_input = ", labor_input)
    run_tests(labor_input)
end
