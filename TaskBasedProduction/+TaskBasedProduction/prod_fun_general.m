function [q, xbar, fval, exitflag] = prod_fun_general(l, z, b_g, e_h)
    % Calculates the quantity produced (q) and task thresholds (xbar)
    % given labor inputs (l), productivity z, and general blueprint density 
    % function and a vector of efficiency functions, one for each labor type.
    %
    % Inputs:
    % - l: Array of labor inputs of different types.
    % - z: Productivity parameter.
    % - b_g: Blueprint density function
    % - e_h: Vector of efficiency functions, one for each type
    %
    % Returns:
    % - q: Quantity produced.
    % - xbar: Array of task thresholds.

    % Objective function for optimization
    function err = objFun(x)
        imp_q = exp(x(1));
        imp_xbar = cumsum(exp(x(2:end)));
        imp_l = imp_q * TaskBasedProduction.unitInputDemand_general(imp_xbar, z, b_g, e_h);
        err = log(imp_l ./ l);
        err = sum(abs(err));  % fminunc requires a single value to minimize
    end

    % Initial guess for optimization
    initial_guess = zeros(length(l), 1);

    % Optimization using fminunc
    options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'off', 'MaxFunctionEvaluations', 1000);
    [x_opt, fval, exitflag] = fminunc(@objFun, initial_guess, options);

    % Check for convergence
    if exitflag <= 0 || max(abs(objFun(x_opt))) > 1e-3
        error('prod_fun: could not find optimal allocation.');
    end

    % Calculate q and xbar
    q = exp(x_opt(1));
    xbar = cumsum(exp(x_opt(2:end)));
end