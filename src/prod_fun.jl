"""
    prod_fun(labor_input::AbstractArray{<:Real}, θ::Real, κ::Real, z::Real, αVec::AbstractArray{<:Real}; initial_guess=nothing, optim_options=nothing)

Calculates the quantity produced (q), and task thresholds (xT) given labor inputs (l), blueprint scale θ, blueprint shape κ, productivity z, and an array of 
comparative advantage values αVec with H elements (one for each worker type).

Inputs:
- `labor_input`: Array of labor inputs of different types.
- `θ`: Blueprint scale parameter.
- `κ`: Blueprint shape parameter.
- `z`: Productivity parameter.
- `αVec`: Array of comparative advantage values with H elements.
- `initial_guess`: (optional) Initial guess for optimization. If not provided, defaults to zeros array.
- `optim_options`: (optional) Optimization options. If not provided, defaults to high tolerance values.

Returns:
- `q`: Quantity produced.
- `xT`: Array of task thresholds.
- `fval`: Final value of the objective function.
"""
function prod_fun(
    labor_input::AbstractArray{<:Real}, 
    θ::Real, 
    κ::Real, 
    z::Real, 
    αVec::AbstractArray{<:Real}; 
    initial_guess=zeros(length(labor_input)), 
    x_tol=1e-8, 
    f_tol=1e-8, 
    g_tol=1e-8, 
    iterations=1000, 
    max_retries=5
)

    function residuals(x::AbstractVector)
        imp_q = exp(x[1])
        imp_xT = cumsum(exp.(x[2:end]))
        imp_l = unitInputDemand(imp_xT, imp_q , θ, κ, z, αVec, true)
        err = log.(imp_l ./ labor_input)
        return err  # Return the error vector for least squares optimization
    end

    retry_count = 0
    while retry_count < max_retries
        try
            # Call the optimizer with Levenberg-Marquardt algorithm and specified options
            result = optimize(residuals, initial_guess, LevenbergMarquardt(), x_tol=x_tol, f_tol=f_tol, g_tol=g_tol, iterations=iterations)
            x_opt = result.minimizer
            fval = maximum(abs.(residuals(x_opt)))

            if fval <= 1e-8
                q = exp(x_opt[1])
                xT = cumsum(exp.(x_opt[2:end]))
                return q, xT, fval
            else
                throw(ErrorException("prod_fun: could not find optimal allocation."))
            end
        catch
            retry_count += 1
            println("prod_fun: An error occurred or could not find optimal allocation. Retrying with a new initial guess. Retry count: $retry_count")
            initial_guess = find_initial_guess(θ, κ, z, αVec; threshold=1e-2)  # Adjust this to your actual method of finding a new initial guess
        end
    end

    error("prod_fun: could not find optimal allocation after $max_retries retries.")
end
