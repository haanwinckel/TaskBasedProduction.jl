"""
prod_fun_general(l::AbstractArray{<:Real}, z::Real, b_g:: Function, e_h::Vector{Function})

Calculates the quantity produced (q), and task thresholds (xbar)
given labor inputs (l),  productivity z, and general blueprint density function and a vector of efficiency functions, one for each labor type.

Inputs:
- l: Array of labor inputs of different types.
- z: Productivity parameter.
- b_g: Blueprint density function
-e_h: Vector of efficiency functions, one for each type 

Returns:
- q: Quantity produced.
- xbar: Array of task thresholds.
"""
function prod_fun_general(l::AbstractArray{<:Real}, z::Real, b_g:: Function, e_h::Vector{Function})

    function objFun(x)
        imp_q = exp(x[1])
        imp_xbar = cumsum(exp.(x[2:end]))
        imp_l = imp_q * unitInputDemand_general(imp_xbar, z,b_g, e_h)
        err = log.(imp_l ./ l)
        return sum(abs.(err))  # Optim requires a single value to minimize
    end

    initial_guess = zeros(length(l))  # Initial guess for optimization
    result = optimize(objFun, initial_guess)
    x_opt = result.minimizer

    if maximum(abs.(objFun(x_opt))) > 1e-4
        error("prod_fun: could not find optimal allocation.")
    end

    q = exp(x_opt[1])
    xbar = cumsum(exp.(x_opt[2:end]))
    return q, xbar
end