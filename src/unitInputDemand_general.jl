
"""
    unitInputDemand_general(xT::Vector{Float64}, z::Real, b_g::Function, e_h::Vector{Function}) -> Vector{Float64}

Calculates unit labor demands given an array `xT` of H-1 thresholds in task space, a productivity value `z`, 
a density function `b_g` for the task distribution, and an array `e_h` of H functions
representing the cost of each labor type as a function of task complexity.

The function first verifies that `b_g` is a valid density function. Then it computes
the labor demand for each labor type by numerically integrating the ratio `b_g(x) / (z * e_h[h](x))`
over the intervals defined by the thresholds in `xT`.

# Arguments
- `xT`: A vector of H-1 thresholds in task space.
- `z`: Productivity value.
- `b_g`: A density function for the task distribution.
- `e_h`: A vector of H functions representing the cost of each labor type as a function of task complexity.

# Returns
- A vector representing the labor demand for each labor type.
"""
function unitInputDemand_general(xT::Vector{Float64}, z::Real, b_g:: Function, e_h::Vector{Function})
    H = length(xT) + 1
    xT = [0.0; xT; Inf]  # Adding 0 and infinity to thresholds
    labor_input = zeros(H)
    # Check if b_g is a density function
    if !is_density_function(b_g, xT[1], xT[end])
        error("b_g(x) is not a valid density function")
    end

    # Compute labor demand for each type
    for h in 1:H
        labor_input[h], _ = quadgk(x -> b_g(x) / (z*e_h[h](x)), xT[h], xT[h+1])
    end

    return labor_input
end



