"""
unitInputDemand(θ::Real, κ::Real, z::Real, αVec::Array{<:Real}, xT::Array{<:Real})

Calculates unit labor demands given blueprint scale θ, blueprint shape κ,
productivity z, an array of comparative advantage values αVec with H elements
(one for each worker type), and an array xT of H-1 thresholds in task space.
"""
function unitInputDemand(θ::Real, κ::Real, z::Real, αVec::AbstractArray{<:Real},
                         xT::AbstractArray{<:Real}, skipParamChecks::Bool = false)
    
    if !skipParamChecks
        @assert θ > 0
        @assert κ > 0
        @assert z > 0
        @assert all(diff(αVec) .> 0) # Alpha values are increasing
        @assert all(diff(xT) .>= 0) # Threshold values are non-decreasing
    end

    H = length(xT) + 1
    l = fill(0.0, H, 1)
    xT = vcat(0.0, xT, Inf) # Add lowest and highest thresholds
    for h in 1:H
        α = αVec[h]
        Υ = α + 1/θ
        if Υ > 0.0
            l[h] = component_positive_ups(Υ, κ, xT[h], xT[h+1]) / Υ^κ
        elseif Υ == 0.0
            l[h] = (xT[h+1]^κ - xT[h]^κ) / (κ * gamma(κ))
        else
            l[h] = component_negative_ups(Υ, κ, xT[h], xT[h+1])
        end
    end
    l = l ./ (z * θ^κ)
end

