"""
    elasticity_sub_comp_general(labor_input::AbstractArray{<:Real}, z::Real, b_g::Function, e_h::Vector{Function}; MPL=nothing, xT=nothing) -> (AbstractArray{<:Real}, AbstractArray{<:Real})

Calculates the elasticity of substitution and complementarity for a given set of parameters.

# Arguments
- `labor_input`: An array of labor inputs of different types with H elements.
- `z`: Productivity parameter.
- `b_g`: General task density function.
- `e_h`: Vector of comparative advantage functions.
- `MPL`: (optional) An array representing the marginal productivity of labor. If not provided, it will be computed within the function.
- `xT`: (optional) An array representing precomputed task thresholds. If not provided, it will be computed within the function.

# Returns
- `ϵ_h_sub`: Matrix of elasticity of substitution values for each worker type h (rows) relative to worker type h_prime (columns).
- `ϵ_h_compl`: Matrix of elasticity of complementarity values for each worker type h (rows) relative to worker type h_prime (columns).
"""
function elasticity_sub_comp_general(labor_input::AbstractArray{<:Real}, z::Real, b_g::Function, e_h::Vector{Function}, MPL=nothing, xT=nothing)
    if xT === nothing
        q, xT = prod_fun_general(labor_input, z, b_g, e_h)
    else
        q, _ = prod_fun_general(labor_input, z, b_g, e_h)  # Recompute q with the given xT
    end

    if MPL === nothing
        MPL = margProdLabor_general(labor_input, z, b_g, e_h)
    end

    H = length(labor_input)
    ρ_h = zeros(Float64, H)
    s_h = MPL .* labor_input / q
    ϵ_h_sub = zeros(Float64, H, H)
    ϵ_h_compl = zeros(Float64, H, H)
    xT = vcat(xT, Inf)  # Add highest thresholds for highest worker type

    e_h_T = [f(xT[i]) for (i, f) in enumerate(e_h[1:end])]
    b_g_T = [b_g(xT[i])/z for i in 1:H]
    
    for h in 1:H-1
        # Compute the numerical derivative of log(e_{h+1} / e_h) with respect to x and evaluate at xT[h]
        log_expr = x -> log(e_h[h+1](x) / e_h[h](x))
        log_derivative = numerical_derivative(log_expr, xT[h])
        
        ρ_h[h] = b_g_T[h] * MPL[h] * (1 / e_h_T[h]) * (1 / log_derivative)
    end

    for h in 1:H
        for h_prime in h+1:H
            if h < h_prime
                if h_prime == h + 1
                    ϵ_h_sub[h, h_prime] = ρ_h[h] / (s_h[h] * s_h[h_prime])
                end
            end
        end
    end

    xi = zeros(H, H, H)
    temp = zeros(H, H, H)

    for h in 1:H
        for h_prime in h+1:H
            for h_bold in 1:H
                xi[h, h_prime, h_bold] = ((h >= h_bold + 1 ? 1 : 0) - sum(s_h[h_bold + 1:H])) * ((h_bold >= h_prime ? 1 : 0) - sum(s_h[1:h_bold]))
                temp[h, h_prime, h_bold] = xi[h, h_prime, h_bold] * (1 / ρ_h[h_bold])
            end
            if h < h_prime
                ϵ_h_compl[h, h_prime] = sum(temp[h, h_prime, 1:H-1])
            else
                ϵ_h_compl[h, h_prime] = NaN
            end
        end
    end

    return ϵ_h_sub, ϵ_h_compl
end
