"""
    margProdLabor_general(labor_input::Union{AbstractArray{<:Real}, Nothing}, z::Real, b_g::Function, e_h::Vector{Function}; xT=nothing, q=nothing) -> AbstractArray{<:Real}

Calculates the marginal productivity of labor for each worker type given the input parameters.

# Arguments
- `labor_input`: An array of labor inputs of different types with H elements. If `nothing`, it will be computed internally given xT and q.
- `z`: A productivity scalar.
- `b_g`: A task density function.
- `e_h`: A vector of comparative advantage functions.
- `xT`: (optional) An array representing the precomputed task thresholds. If not provided, it will be computed within the function.
- `q`: (optional) A scalar representing the precomputed quantity produced. If not provided, it will be computed within the function.

# Returns
- An array representing the marginal productivity of labor for each worker type.

If `labor_input` is not provided, it will be computed using the `q` and `unitInputDemand_general` function based on the other parameters.
"""
function margProdLabor_general(
    labor_input::Union{AbstractArray{<:Real}, Nothing},
    z::Real,
    b_g::Function,
    e_h::Vector{Function},
    xT=nothing,
    q=nothing)
    # Calculate q and xT if they are not provided
    if xT === nothing || q === nothing
        q, xT = prod_fun_general(labor_input, z, b_g, e_h)
    end
    
    # If labor_input is missing, calculate it
    if labor_input === nothing
        labor_input = q * unitInputDemand_general(xT, q, z, b_g, e_h)
    end

    H = length(e_h)
    mpl_over_mpl1 = [1.0]
    
    # Calculate the ratio e_{h} / e_{h-1} for h = 2:H and evaluate at xT[h-1]
    temp = zeros(H-1)
    for h in 2:H
        ratio_value = e_h[h](xT[h-1]) / e_h[h-1](xT[h-1])
        temp[h-1] = ratio_value
    end
    mpl_over_mpl1 = [1; cumprod(temp, dims=1)]
    mpl1 = q / sum(mpl_over_mpl1 .* labor_input)
    mpl = mpl_over_mpl1 * mpl1

    return mpl
end

