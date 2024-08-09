"""
    margProdLabor(labor_input::Union{AbstractArray{<:Real}, Nothing}, θ::Real, κ::Real, z::Real, αVec::AbstractArray{<:Real}; xT=nothing, q=nothing) -> AbstractArray{<:Real}

Calculates the marginal productivity of labor for each worker type given the input parameters.

# Arguments
- `labor_input`: An array of labor demand values. If `nothing`, it will be computed internally (given xT and q).
- `θ`: Blueprint scale parameter.
- `κ`: Blueprint shape parameter.
- `z`: Productivity parameter.
- `αVec`: An array of comparative advantage values.
- `xT`: (optional) An array representing the precomputed task thresholds. If not provided, it will be computed within the function.
- `q`: (optional) A scalar representing the precomputed quantity produced. If not provided, it will be computed within the function.

# Returns
- An array representing the marginal productivity of labor for each worker type.

If `labor_input` is not provided, it will be computed using the `q` and `unitInputDemand` function based on the other parameters.
"""
function margProdLabor(
    labor_input::Union{AbstractArray{<:Real}, Nothing},
    θ::Real,
    κ::Real,
    z::Real,
    αVec::AbstractArray{<:Real},
    xT=nothing,
    q=nothing)
    # Calculate q and xT if they are not provided
    if xT === nothing || q === nothing
        q, xT = prod_fun(labor_input, θ, κ, z, αVec)
    end
    
    # If labor_input is missing, calculate it using unitInputDemand
    if labor_input === nothing
        labor_input = unitInputDemand(xT, q, θ, κ, z, αVec)
    end

    # Compute marginal productivity of labor
    mpl_over_mpl1 = [1; cumprod(exp.(diff(αVec, dims=1) .* xT), dims=1)]
    mpl1 = q / sum(mpl_over_mpl1 .* labor_input)
    mpl = mpl_over_mpl1 * mpl1

    return mpl
end


