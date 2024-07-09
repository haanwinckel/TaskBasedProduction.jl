function val = component_positive_ups(upsilon, kappa, xT_low, xT_high)
    % component_positive_ups Computes the positive ups component
    %
    % Inputs:
    %   upsilon - parameter upsilon (scalar)
    %   kappa - parameter kappa (scalar)
    %   xT_low - lower threshold (scalar)
    %   xT_high - upper threshold (scalar)
    %
    % Output:
    %   val - computed value (scalar)
    
    if isinf(xT_high)
        upper_element = 1;
    else
        upper_element = gammainc(upsilon * xT_high, kappa, 'lower');
    end

    if xT_low == 0.0
        lower_element = 0.0;
    else
        lower_element = gammainc(upsilon * xT_low, kappa, 'lower');
    end
    val = upper_element - lower_element;
end