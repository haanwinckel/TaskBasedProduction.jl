"""
    elasticity_sub_comp(labor_input::AbstractArray{<:Real}, θ::Real, κ::Real, z::Real, αVec::AbstractArray{<:Real}) -> (AbstractArray{<:Real}, AbstractArray{<:Real})

Calculates the elasticity of substitution and complementarity for a given set of parameters.

# Arguments
- `labor_input`: An array of labor inputs of different types with H elements.
- `θ`: Blueprint scale parameter.
- `κ`: Blueprint shape parameter.
- `z`: Productivity parameter.
- `αVec`: An array of comparative advantage values with H elements.

# Returns
- `ϵ_h_sub`: Matrix of elasticity of substitution values for each worker type h (rows) relative to worker type h_prime (columns).
- `ϵ_h_compl`: Matrix of elasticity of complementarity values for each worker type h (rows) relative to worker type h_prime (columns).
"""
function elasticity_sub_comp(labor_input::AbstractArray{<:Real}, θ::Real, κ::Real, z::Real, αVec::AbstractArray{<:Real})
    q, xT= prod_fun(labor_input, θ, κ, z, αVec)
    MPL= margProdLabor(labor_input, θ, κ, z, αVec)
    H = length(αVec)
    ρ_h = zeros(Float64, H)
    s_h = MPL.*labor_input/q
    ϵ_h_sub=zeros(Float64, H, H)
    ϵ_h_compl=zeros(Float64, H, H)
    xT = vcat(xT, Inf)  # Add highest thresholds for highest worker type
    e_h_T = exp.(αVec .* xT)  # Denominator of ρ_h



    for h in 1:H-1
        b_g_T = xT[h]^(κ - 1) * (1 / (z* gamma(κ) * θ^κ)) * exp(-xT[h] / θ)
        ρ_h[h] = b_g_T * MPL[h] * (1 / e_h_T[h]) * (1 / (αVec[h+1] - αVec[h]))
    end

    for h in 1:H
        for h_prime in h+1:H
            if h<h_prime
                if h_prime == h + 1
                ϵ_h_sub[h, h_prime] = ρ_h[h] / (s_h[h] * s_h[h_prime])
                end
            end

        end
    end
    xi = zeros(H, H, H)
    temp= zeros(H,H, H)

    for h in 1:H
        for h_prime in h+1:H
            for h_bold in 1:H 
                xi[h, h_prime, h_bold]=((h>=h_bold+1 ? 1 : 0)-sum(s_h[h_bold+1:H]))*((h_bold>=h_prime ? 1 : 0)-sum(s_h[1:h_bold]))
                temp[h, h_prime, h_bold]=xi[h, h_prime, h_bold]*(1/ρ_h[h_bold])
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
