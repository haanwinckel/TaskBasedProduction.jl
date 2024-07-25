"""
    find_initial_guess_gen(labor_input::AbstractArray{<:Real}, z::Real, αVec::AbstractArray{<:Real}, pdf::Function; threshold::Real=1e-2, verbose::Bool=false)

Generate an initial guess for the optimization problem using a general density function such that the implied labor demand is non-trivial.

# Arguments
- `labor_input::AbstractArray{<:Real}`: The observed labor input for each task.
- `z::Real`: A scaling factor for the labor input.
- `αVec::AbstractArray{<:Real}`: An array of task-specific parameters.
- `pdf::Function`: The general density function.
- `threshold::Real`: The minimum acceptable labor demand for each task.
- `verbose::Bool`: Optional boolean flag to enable or disable verbose output for debugging.

# Returns
- `initial_guess::Array{<:Real}`: A vector containing the initial guess for the optimization, including the log of the initial production quantity `q` and the initial task thresholds `xT`.

# Description
This function generates an initial guess for the optimization problem by:
1. Fixing the initial guess for `q` at 1.
2. Generating initial `xT` values using random percentiles from the provided CDF function.
3. Adjusting the `xT` values iteratively to ensure the implied labor demand for each task is above the specified threshold.

If the implied labor demand for any task is below the threshold, the `xT` values are re-shuffled using the `generate_initial_xT` function. This process continues until the implied labor demand for all tasks is above the threshold or the maximum number of iterations is reached.

If the adjustment process encounters an error, new `xT` values are generated from scratch.

Verbose output can be enabled by setting the `verbose` parameter to `true`, which will print debug information during the percentile calculation.
"""
function find_initial_guess_gen(z::Real, b_g::Function, e_h::Vector{Function}; threshold::Real=1e-2, verbose::Bool=false)
    H = length(e_h)  # Number of tasks
    # Initial guess for q is fixed at 1
    initial_q = 0.0  # log(1) is 0
    
    # Define the CDF by integrating the PDF
    function c_cdf(b_g::Function, x::Real)
        integral, _ = quadgk(b_g, 0, x)
        return integral
    end

    # Function to calculate the percentile using the general CDF
    function calculate_percentile(p::Real, b_g::Function)
        # Define the error function for optimization
        function err(x::Vector{Float64})
            if verbose
                println("Evaluating CDF at x[1] = $(x[1])")
            end
            if x[1] < 0
                return Inf  # Return a large error for invalid domain values
            end
            cdf_val = c_cdf(b_g, x[1])
            if verbose
                println("CDF value: $cdf_val, Target percentile: $p")
            end
            return (cdf_val - p)^2  # Squared difference to ensure non-negative optimization
        end

        # Set the initial guess for optimization
        initial_guess = [0.1]  # Initial guess as a vector
        # Perform the optimization using the BFGS method
        result = optimize(err, initial_guess, BFGS())
        # Extract the percentile value
        x_percentile = Optim.minimizer(result)[1]

        return x_percentile
    end

    # Function to generate initial xT values using random percentiles from the general CDFs
    function generate_initial_xT()
        percentiles = sort(rand(H-1))  # Generate random percentiles and sort them
        xT = [calculate_percentile(p, b_g) for p in percentiles]  # Calculate xT for each percentile
        return xT
    end

    function adjust_xT(xT)
        imp_l = zeros(H)  # Placeholder for labor demand
        max_iterations = 1000
        iteration = 0

        while any(imp_l .< threshold) && iteration < max_iterations
            try
                imp_xT = cumsum(exp.(xT[1:end]))
                imp_l = exp(initial_q) * unitInputDemand_general(imp_xT, z, b_g, e_h)
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
            error("find_initial_guess_gen: Could not find suitable initial xT within the maximum iterations. Consider changing tolerance.")
        end

        return xT
    end
    
    initial_xT = generate_initial_xT()
    adjusted_xT = adjust_xT(initial_xT)

    # Combine q and xT into the initial guess vector
    initial_guess = vcat([initial_q], adjusted_xT)
    
    return initial_guess
end