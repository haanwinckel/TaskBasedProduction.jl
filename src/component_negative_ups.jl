using SpecialFunctions: gamma
function component_negative_ups(Υ::Real, κ::Real, xT_low::Real, xT_high::Real)
    MAX_ITER = 100
    TOL = 1e-16

    m = 0
    upper_element = xT_high^κ * exp(-Υ * xT_high)
    lower_element = xT_low^κ * exp(-Υ * xT_low)
    val = (upper_element - lower_element) / gamma(κ + 1)

    converged = false
    change = 0
    while !converged && m < MAX_ITER
        m += 1
        upper_element *= Υ * xT_high
        lower_element *= Υ * xT_low
        change = (upper_element - lower_element) / gamma(κ + m + 1)

        m += 1
        upper_element *= Υ * xT_high
        lower_element *= Υ * xT_low
        change += (upper_element - lower_element) / gamma(κ + m + 1)

        converged = abs(change) < TOL
        val += change
    end
    if !converged
        return NaN
    else
        return val
    end
end