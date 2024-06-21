"""
margProdLabor(inputDemand::Array{<:Real}, αVec::Array{<:Real}, xT::Array{<:Real})

Calculate marginal products of labor for each worker type given an array
of H labor demand values, an array of comparative advantage values α, and
an array of H-1 task thresholds xT that corresponds to that labor demand.
"""
function margProdLabor(inputDemand::AbstractArray{<:Real}, αVec::AbstractArray{<:Real}, xT::AbstractArray{<:Real})
    mpl_over_mpl1 = [1; cumprod(exp.(diff(αVec, dims=1) .* xT), dims=1)]
    mpl1 = 1 / sum(mpl_over_mpl1 .* inputDemand)
    mpl = mpl_over_mpl1 * mpl1
    return mpl
end

"""
margProdLabor(θ::Real, κ::Real, z::Real, αVec::Array{<:Real}, xT::Array{<:Real})

Calculate marginal products of labor for each worker type given blueprint
characteristics (scale θ, shape κ, productivity z), an array of comparative 
advantage values α, and an array of H-1 task thresholds xT.

NOTE: if you already have calculated the labor demand corresponding to
these inputs, it is more efficient to use the alternative version of this
function that takes labor demands as inputs.
"""
function margProdLabor(θ::Real, κ::Real, z::Real, αVec::Array{<:Real}, xT::Array{<:Real})
    l = unitInputDemand(θ, κ, z, αVec, xT)
    return margProdLabor(l, αVec, xT)
end