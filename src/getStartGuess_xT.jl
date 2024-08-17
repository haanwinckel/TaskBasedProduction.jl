"""
    getStartGuess_xT(labor_input::AbstractArray{<:Real}, θ::Real, κ::Real, z::Real, αVec::AbstractArray{<:Real}; threshold::Real=1e-2)

Generate an initial guess for the optimization problem in `prod_fun` such that the implied labor demand is non-trivial.

# Arguments
- `labor_input::AbstractArray{<:Real}`: The observed labor input for each task.
- `θ::Real`: The scale parameter of the gamma distribution.
- `κ::Real`: The shape parameter of the gamma distribution.
- `z::Real`: A scaling factor for the labor input.
- `αVec::AbstractArray{<:Real}`: An array of task-specific parameters.
- `threshold::Real`: The minimum acceptable labor demand for each task.

# Returns
- `initial_guess::Array{<:Real}`: A vector containing the initial guess for the optimization, including the log of the initial production quantity `q` and the initial task thresholds `xT`.

# Description
This function generates an initial guess for the optimization problem by:
1. Fixing the initial guess for `q` at 1.
2. Generating initial `xT` values using random percentiles from the gamma distribution defined by `θ` and `κ`.
3. Adjusting the `xT` values iteratively to ensure the implied labor demand for each task is above the specified threshold.

If the implied labor demand for any task is below the threshold, the `xT` values are re-shuffled using the `generate_initial_xT` function. This process continues until the implied labor demand for all tasks is above the threshold or the maximum number of iterations is reached.

If the adjustment process encounters an error, new `xT` values are generated from scratch.
"""
function  getStartGuess_xT(θ::Real, κ::Real, z::Real, αVec::AbstractArray{<:Real}; threshold::Real=1e-2)
    H = length(αVec)  # Number of tasks
    # Initial guess for q is fixed at 1
    initial_q = 0.0  # log(1) is 0
    
    function generate_initial_xT()
        # Generate initial xT values using random percentiles from the gamma distribution
        gamma_dist = Gamma(κ, θ)
        percentiles = sort(rand(H-1))  # Generate random percentiles and sort them
        xT = [quantile(gamma_dist, p) for p in percentiles]  # Calculate xT for each percentile
        return xT
    end
    
    function adjust_xT(xT)
        imp_l = zeros(H)  # Placeholder for labor demand
        max_iterations = 1000
        iteration = 0

        while any(imp_l .< threshold) && iteration < max_iterations
            try
                imp_xT = cumsum(exp.(xT[1:end]))
                imp_q=exp(initial_q)
                imp_l = unitInputDemand(imp_xT, imp_q, θ, κ, z, αVec)
            catch
                # If there's an error, generate new initial xT values from scratch
                xT = generate_initial_xT()
                continue
            end

            if any(imp_l .< threshold)
                xT = generate_initial_xT()
            end

            iteration += 1
        end

        if iteration == max_iterations
            error("find_initial_guess: Could not find suitable initial xT within the maximum iterations. Consider changing tolerance.")
        end

        return xT
    end
    
    initial_xT = generate_initial_xT()
    adjusted_xT = adjust_xT(initial_xT)

    # Combine q and xT into the initial guess vector
    initial_guess = vcat([initial_q], adjusted_xT)
    
    return initial_guess
end
