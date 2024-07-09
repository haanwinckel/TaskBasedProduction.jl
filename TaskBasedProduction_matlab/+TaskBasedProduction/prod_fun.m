function [q, xT] = prod_fun(l, theta, kappa, z, alphaVec)
    % prod_fun Calculates the quantity produced and task thresholds
    %
    % Calculates the quantity produced (q), and task thresholds (xT)
    % given labor inputs (l), blueprint scale theta, blueprint shape kappa,
    % productivity z, and an array of comparative advantage values alphaVec
    % with H elements (one for each worker type).
    %
    % Inputs:
    %   l - array of labor inputs of different types (vector)
    %   theta - blueprint scale parameter (scalar)
    %   kappa - blueprint shape parameter (scalar)
    %   z - productivity parameter (scalar)
    %   alphaVec - array of comparative advantage values with H elements (vector)
    %
    % Returns:
    %   q - quantity produced (scalar)
    %   xT - array of task thresholds (vector)
    
    % Objective function for optimization
    function val = objFun(x)
        imp_q = exp(x(1));
        imp_xT = cumsum(exp(x(2:end)));
        imp_l = imp_q * TaskBasedProduction.unitInputDemand(theta, kappa, z, alphaVec, imp_xT, true);
        err = log(imp_l ./ l);
        val = sum(abs(err));  % Optimization requires a single value to minimize
    end

    % Initial guess for optimization
    initial_guess = zeros(length(alphaVec), 1);
    
    % Set options for fminunc
    options = optimoptions('fminunc', 'Display', 'off', 'Algorithm', 'quasi-newton');
    
    % Perform optimization
    [x_opt, ~, exitflag] = fminunc(@objFun, initial_guess, options);

    % Check if optimization was successful
    if exitflag <= 0
        error('prod_fun: could not find optimal allocation.');
    end

    % Calculate q and xT
    q = exp(x_opt(1));
    xT = cumsum(exp(x_opt(2:end)));
end