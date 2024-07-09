function valid = is_density_function(b_g, a, b)
    % Check if b_g integrates to 1 over [a, b]
    integral_value = integral(b_g, a, b);
    valid = abs(integral_value - 1) < 1e-6;
end