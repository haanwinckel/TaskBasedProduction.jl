"""
prod_fun(l::AbstractArray{<:Real}, θ::Real, κ::Real, z::Real, αVec::AbstractArray{<:Real})

Calculates the quantity produced (q), and task thresholds (xT)
given labor inputs (l), blueprint scale θ, blueprint shape κ, productivity z, and an array of 
comparative advantage values αVec with H elements (one for each worker type).

Inputs:
- l: Array of labor inputs of different types.
- θ: Blueprint scale parameter.
- κ: Blueprint shape parameter.
- z: Productivity parameter.
- αVec: Array of comparative advantage values with H elements.

Returns:
- q: Quantity produced.
- xT: Array of task thresholds.
"""
function prod_fun(labor_input::AbstractArray{<:Real}, θ::Real, κ::Real, z::Real, αVec::AbstractArray{<:Real})

    function objFun(x)
        imp_q = exp(x[1])
        imp_xT = cumsum(exp.(x[2:end]))
        imp_l = imp_q * unitInputDemand( imp_xT, θ, κ, z, αVec, true)
        err = log.(imp_l ./ labor_input)
        return sum(abs.(err))  # Optim requires a single value to minimize
    end

    initial_guess = zeros(length(αVec))  # Initial guess for optimization
    result = optimize(objFun, initial_guess, Optim.Options(g_tol=1e-12, f_tol=1e-12, x_tol=1e-12, iterations=1000))
    x_opt = result.minimizer
   if maximum(abs.(objFun(x_opt))) > 1e-4
       error("prod_fun: could not find optimal allocation.")
   end
    fval=maximum(abs.(objFun(x_opt)))
    q = exp(x_opt[1])
    xT = cumsum(exp.(x_opt[2:end]))
    return q, xT, fval
end

