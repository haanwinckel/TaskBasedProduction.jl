function val = component_negative_ups(upsilon, kappa, xT_low, xT_high)
    % component_negative_ups Computes the negative ups component
    %
    % Inputs:
    %   upsilon - parameter upsilon (scalar)
    %   kappa - parameter kappa (scalar)
    %   xT_low - lower threshold (scalar)
    %   xT_high - upper threshold (scalar)
    %
    % Output:
    %   val - computed value (scalar)
    
    MAX_ITER = 100;
    TOL = 1e-16;

    m = 0;
    upper_element = xT_high^kappa * exp(-upsilon * xT_high);
    lower_element = xT_low^kappa * exp(-upsilon * xT_low);
    val = (upper_element - lower_element) / gamma(kappa + 1);

    converged = false;
    change = 0;
    while ~converged && m < MAX_ITER
        m = m + 1;
        upper_element = upper_element * upsilon * xT_high;
        lower_element = lower_element * upsilon * xT_low;
        change = (upper_element - lower_element) / gamma(kappa + m + 1);

        m = m + 1;
        upper_element = upper_element * upsilon * xT_high;
        lower_element = lower_element * upsilon * xT_low;
        change = change + (upper_element - lower_element) / gamma(kappa + m + 1);

        converged = abs(change) < TOL;
        val = val + change;
    end
    if ~converged
        val = NaN;
    end
end