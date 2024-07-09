function l = unitInputDemand(theta, kappa, z, alphaVec, xT, skipParamChecks)
    % unitInputDemand Calculates unit labor demands
    % 
    % Calculates unit labor demands given blueprint scale theta, blueprint shape kappa,
    % productivity z, an array of comparative advantage values alphaVec with H elements
    % (one for each worker type), and an array xT of H-1 thresholds in task space.
    %
    % Inputs:
    %   theta - blueprint scale (scalar)
    %   kappa - blueprint shape (scalar)
    %   z - productivity (scalar)
    %   alphaVec - array of comparative advantage values (vector)
    %   xT - array of thresholds in task space (vector)
    %   skipParamChecks - flag to skip parameter checks (boolean)
    %
    % Output:
    %   l - unit labor demands (vector)

    % Assign default value to skipParamChecks if not provided
    if nargin < 6
        skipParamChecks = false;
    end

    if ~skipParamChecks
        assert(theta > 0, 'theta must be greater than 0');
        assert(kappa > 0, 'kappa must be greater than 0');
        assert(z > 0, 'z must be greater than 0');
        assert(all(diff(alphaVec) > 0), 'alphaVec values must be increasing');
        assert(all(diff(xT) >= 0), 'xT values must be non-decreasing');
    end

    H = length(xT) + 1;
    l = zeros(H, 1);
    xT = [0.0; xT(:); Inf]; % Add lowest and highest thresholds

    for h = 1:H
        alpha = alphaVec(h);
        upsilon = alpha + 1/theta;
        if upsilon > 0.0
            l(h) = TaskBasedProduction.component_positive_ups(upsilon, kappa, xT(h), xT(h+1)) / upsilon^kappa;
        elseif upsilon == 0.0
            l(h) = (xT(h+1)^kappa - xT(h)^kappa) / (kappa * gamma(kappa));
        else
            l(h) = TaskBasedProduction.component_negative_ups(upsilon, kappa, xT(h), xT(h+1));
        end
    end
    l = l / (z * theta^kappa);
end
