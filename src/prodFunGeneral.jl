"""
    prodFunGeneral(labor_input::AbstractArray{<:Real}, z::Real, b_g:: Function, e_h::Vector{Function}; initial_guess=nothing, x_tol=1e-12, f_tol=1e-12, g_tol=1e-12, iterations=1000, max_retries=5)

Calculates the quantity produced (q), and task thresholds (xT) given labor inputs (labor_input), productivity z, general blueprint density function (b_g), and a vector of efficiency functions (e_h), one for each labor type.

Inputs:
- `labor_input`: Array of labor inputs of different types.
- `z`: Productivity parameter.
- `b_g`: Blueprint density function.
- `e_h`: Vector of efficiency functions, one for each type.
- `initial_guess`: (optional) Initial guess for optimization. If not provided, defaults to zeros array.
- `x_tol`: (optional) Tolerance for the solution vector. Default is 1e-12.
- `f_tol`: (optional) Tolerance for the function value. Default is 1e-12.
- `g_tol`: (optional) Tolerance for the gradient. Default is 1e-12.
- `iterations`: (optional) Maximum number of iterations for the optimization. Default is 1000.
- `max_retries`: (optional) Maximum number of retries if the optimization fails. Default is 5.

Returns:
- `q`: Quantity produced.
- `xT`: Array of task thresholds.
"""
function prodFunGeneral(
    labor_input::AbstractArray{<:Real}, 
    z::Real, 
    b_g::Function, 
    e_h::Vector{Function}; 
    initial_guess=nothing,  # Allow for a default of nothing
    x_tol=1e-12, 
    f_tol=1e-12, 
    g_tol=1e-12, 
    iterations=1000, 
    max_retries=5
)

    # If initial_guess is nothing, calculate it using getStartGuessGen_xT
    if initial_guess === nothing
        initial_guess = getStartGuessGen_xT(z, b_g, e_h)
    end
    function objFun(x)
        imp_q = exp(x[1])
        imp_xT = cumsum(exp.(x[2:end]))
        imp_l = unitInputDemandGeneral(imp_xT, imp_q, z, b_g, e_h)
        err = log.(imp_l ./ labor_input)
        return err  # Return the error vector for least squares optimization
    end

    retry_count = 0
    while retry_count < max_retries
        try
            # Call the optimizer with Levenberg-Marquardt algorithm and specified options
            result = optimize(objFun, initial_guess, LevenbergMarquardt(), x_tol=x_tol, f_tol=f_tol, g_tol=g_tol, iterations=iterations)
            x_opt = result.minimizer
            fval = maximum(abs.(objFun(x_opt)))

            if fval <= f_tol
                q = exp(x_opt[1])
                xT = cumsum(exp.(x_opt[2:end]))
                return q, xT
            else
                throw(ErrorException("prod_fun_general: could not find optimal allocation."))
            end
        catch
            retry_count += 1
            println("prod_fun_general: An error occurred or could not find optimal allocation. Retrying with a new initial guess. Retry count: $retry_count")
            initial_guess = getStartGuessGen_xT(z, b_g, e_h)  # Adjust this to your actual method of finding a new initial guess
        end
    end

    error("prod_fun_general: could not find optimal allocation after $max_retries retries.")
end
