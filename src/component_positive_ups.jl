function component_positive_ups(Υ::Real, κ::Real, xT_low::Real, xT_high::Real)
    if isinf(xT_high)
        upper_element = 1
    else
        upper_element = low_norm_gamma_inc(κ, Υ * xT_high)
    end

    if xT_low == 0.0
        lower_element = 0.0
    else
        lower_element = low_norm_gamma_inc(κ, Υ * xT_low)
    end
    return upper_element - lower_element
end