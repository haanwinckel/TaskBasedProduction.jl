function numerical_derivative(f::Function, x::Real; δ::Real=1e-5)
    return (f(x + δ) - f(x - δ)) / (2 * δ) #Central difference method to ensure a better approximation
end