function l_h = unitInputDemand_general(xT, z, b_g, e_h)
    % Calculates unit labor demands given an array xT of H-1 thresholds in task space,
    % productivity value z, a density function b_g for the task distribution,
    % and an array e_h of H functions representing the cost of each labor type
    % as a function of task complexity.
    
    % Ensure xT is a column vector
    xT = xT(:);
    
    % Number of labor types
    H = length(xT) + 1;
    
    % Adding 0 and infinity to thresholds
    xT = [0; xT; Inf];
    
    % Initialize the output vector
    l_h = zeros(H, 1);
    
    % Check if b_g is a density function
    if ~TaskBasedProduction.is_density_function(b_g, xT(1), xT(end))
        error('b_g(x) is not a valid density function');
    end
    
    % Compute labor demand for each type
    for h = 1:H
        integrand = @(x) b_g(x) / (z * e_h{h}(x));
        l_h(h) = integral(integrand, xT(h), xT(h+1),  'ArrayValued', true);
    end
end


