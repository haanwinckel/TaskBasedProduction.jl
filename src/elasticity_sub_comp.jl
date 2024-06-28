"""
elasticity_sub_comp(xT::AbstractArray{<:Real}, l::AbstractArray{<:Real}, q::Real, MPL::AbstractArray{<:Real}, θ::Union{Real, AbstractArray{<:Real}}, κ::Union{Real, AbstractArray{<:Real}}, z::Union{Real, AbstractArray{<:Real}}, αVec::AbstractArray{<:Real})

Calculates the elasticity of substitution and complementarity for a given set of parameters.

Parameters:
- xT: Array of task thresholds with H-1 elements.
- l: Array of labor inputs of different types with H elements.
- q: Total output.
- MPL: Array of marginal products of labor for each worker type with H elements.
- θ: Blueprint scale parameter (scalar or array indexed by h with H elements).
- κ: Blueprint shape parameter (scalar or array indexed by h with H elements).
- z: Productivity parameter (scalar or array indexed by h with H elements).
- αVec: Array of comparative advantage values with H elements.

Returns:
- ϵ_h_sub: Matrix of elasticity of substitution values for each worker type h (rows) relative to worker type h_prime (columns).
- ϵ_h_comp: Matrix of elasticity of complementarity values for each worker type h (rows) relative to worker type h_prime (columns).

"""
function elasticity_sub_comp(xT::AbstractArray{<:Real}, l::AbstractArray{<:Real}, q::Real, MPL::AbstractArray{<:Real}, θ::Union{Real, AbstractArray{<:Real}}, κ::Union{Real, AbstractArray{<:Real}}, z::Union{Real, AbstractArray{<:Real}}, αVec::AbstractArray{<:Real})
    H = length(αVec)
    ρ_h = zeros(Float64, H)
    s_h = MPL.*l/q
    ϵ_h_sub=zeros(Float64, H, H)
    ϵ_h_compl=zeros(Float64, H, H)
    xT = vcat(xT, Inf)  # Add highest thresholds for highest worker type
    e_h_T = exp.(αVec .* xT)  # Denominator of ρ_h

    # Helper function to handle scalar or array inputs for θ, κ, and z
    function get_param(param, h)
        return isa(param, AbstractArray) ? param[h] : param
    end

    for h in 1:H-1
        θ_h = get_param(θ, h)
        κ_h = get_param(κ, h)
        z_h = get_param(z, h)
        b_g_T = xT[h]^(κ_h - 1) * (1 / (z_h * gamma(κ_h) * θ_h^κ_h)) * exp(-xT[h] / θ_h)
        ρ_h[h] = b_g_T * MPL[h] * (1 / e_h_T[h]) * (1 / (αVec[h+1] - αVec[h]))
    end

    for h in 1:H
        for h_prime in 1:H
            if h_prime == h + 1
                ϵ_h_sub[h, h_prime] = ρ_h[h] / (s_h[h] * s_h[h_prime])
            end
        end
    end
    xi = zeros(H, H, H)
    temp= zeros(H,H, H)

    for h in 1:H
        for h_prime in 1:H
            for h_bold in 1:H 
                xi[h, h_prime, h_bold]=((h>=h_bold+1 ? 1 : 0)-sum(s_h[h_bold+1:H]))*((h_bold>=h_prime ? 1 : 0)-sum(s_h[1:h_bold]))
                temp[h, h_prime, h_bold]=xi[h, h_prime, h_bold]*(1/ρ_h[h_bold])
            end
            ϵ_h_compl[h, h_prime]=sum(temp[h, h_prime, 1:H-1])
        end
    end

   return ϵ_h_sub, ϵ_h_compl
end
