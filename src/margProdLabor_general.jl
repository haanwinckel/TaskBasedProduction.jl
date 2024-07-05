"""
margProdLabor_general(inputDemand::Array{<:Real}, αVec::Array{<:Real}, xT::Array{<:Real})

Calculate marginal products of labor for each worker type given an array
of H labor demand values, a vector of comparative advantage functions e_h, and
an array of H-1 task thresholds xT that corresponds to that labor demand.
"""
function margProdLabor_general(xT::AbstractArray{<:Real},l::AbstractArray{<:Real}, e_h::Vector{Function})
    H = length(e_h)
    mpl_over_mpl1 = [1.0]
    # Calculate the ratio e_{h} / e_{h-1} for h = 2:H and evaluate at xT[h-1]
    temp = zeros(H-1)
    for h in 2:H
        ratio_value = e_h[h](xT[h-1]) / e_h[h-1](xT[h-1])
        temp[h-1]=ratio_value
    end
    mpl_over_mpl1=[1; cumprod(temp,dims=1)]
    mpl1 = 1 / sum(mpl_over_mpl1 .* l)
    mpl = mpl_over_mpl1 * mpl1
    return mpl
end