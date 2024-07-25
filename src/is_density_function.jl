# Function to check if b_g(x) is a density function
function is_density_function(b_g :: Function, lower::Float64, upper::Float64)
    integral, _ = quadgk(b_g, lower, upper)
    return abs(integral - 1.0) < 1e-6
end