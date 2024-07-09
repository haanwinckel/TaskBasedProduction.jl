"""
elasticity_sub_comp_general(xT::AbstractArray{<:Real},l::AbstractArray{<:Real}, q::Real, MPL::AbstractArray{<:Real}, z:: Real, b_g::Function, e_h::Vector{Function})

Calculates the elasticity of substitution and complementarity for a given set of parameters.

Parameters:
- xT: Array of task thresholds with H-1 elements.
- l: Array of labor inputs of different types with H elements.
- q: Total output.
- MPL: Array of marginal products of labor for each worker type with H elements.
- b_g: general task density function
- e_h: vector of comparative advantage functions

Returns:
- ϵ_h_sub: Matrix of elasticity of substitution values for each worker type h (rows) relative to worker type h_prime (columns).
- ϵ_h_comp: Matrix of elasticity of complementarity values for each worker type h (rows) relative to worker type h_prime (columns).

"""


function elasticity_sub_comp_general(xT::AbstractArray{<:Real},l::AbstractArray{<:Real}, q::Real, MPL::AbstractArray{<:Real}, z:: Real, b_g::Function, e_h::Vector{Function})
    H = length(l)
    ρ_h = zeros(Float64, H)
    s_h = MPL .* l / q
    ϵ_h_sub = zeros(Float64, H, H)
    ϵ_h_compl = zeros(Float64, H, H)
    xT = vcat(xT, Inf)  # Add highest thresholds for highest worker type

    e_h_T = [f(xT[i]) for (i, f) in enumerate(e_h[1:end])]
    b_g_T = [b_g(xT[i])/z for i in 1:H]
    x = symbols("x")
    for h in 1:H-1
                # Compute the derivative of log(e_{h+1} / e_h) with respect to x and evaluate at xT[h]
                log_expr = log(e_h[h+1](x) / e_h[h](x))
                diff_log_expr = diff(log_expr, x)
                log_derivative = diff_log_expr(x => xT[h])
                
                ρ_h[h] = b_g_T[h] * MPL[h] * (1 / e_h_T[h]) * (1/log_derivative)
    end

    for h in 1:H
        for h_prime in 1:H
            if h<h_prime
                if h_prime == h + 1
                ϵ_h_sub[h, h_prime] = ρ_h[h] / (s_h[h] * s_h[h_prime])
                end
            else ϵ_h_sub[h, h_prime] = NaN
            end

        end
    end

    xi = zeros(H, H, H)
    temp = zeros(H, H, H)

    for h in 1:H
        for h_prime in h:H
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