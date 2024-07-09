function [epsilon_h_sub, epsilon_h_compl] = elasticity_sub_comp_general(xT, l, q, MPL, z, b_g, e_h)
    % Calculates the elasticity of substitution and complementarity for a given set of parameters.
    % Inputs:
    %   xT - Array of task thresholds with H-1 elements
    %   l - Array of labor inputs of different types with H elements
    %   q - Total output
    %   MPL - Array of marginal products of labor for each worker type with H elements
    %   z - Some scalar parameter
    %   b_g - General task density function (handle)
    %   e_h - Cell array of comparative advantage functions (cell array of function handles)
    % Outputs:
    %   epsilon_h_sub - Matrix of elasticity of substitution values (H x H)
    %   epsilon_h_compl - Matrix of elasticity of complementarity values (H x H)

    H = length(l);
    rho_h = zeros(H, 1);
    s_h = (MPL .* l) / q;
    epsilon_h_sub = NaN(H, H);
    epsilon_h_compl = NaN(H, H);
    xT = [xT; Inf];  % Add highest threshold for highest worker type

    e_h_T = zeros(H, 1);
    b_g_T = zeros(H, 1);

    for i = 1:H
        e_h_T(i) = e_h{i}(xT(i));
        b_g_T(i) = b_g(xT(i)) / z;
    end

    syms x;
    for h = 1:H-1
        % Compute the derivative of log(e_{h+1} / e_h) with respect to x and evaluate at xT[h]
        log_expr = log(e_h{h+1}(x) / e_h{h}(x));
        diff_log_expr = diff(log_expr, x);
        log_derivative = double(subs(diff_log_expr, x, xT(h)));

        rho_h(h) = b_g_T(h) * MPL(h) * (1 / e_h_T(h)) * (1 / log_derivative);
    end

    for h = 1:H
        for h_prime = 1:H
            if h < h_prime && h_prime == h + 1
                epsilon_h_sub(h, h_prime) = rho_h(h) / (s_h(h) * s_h(h_prime));
            end
        end
    end

    xi = zeros(H, H, H);
    temp = zeros(H, H, H);

    for h = 1:H
        for h_prime = h:H
            for h_bold = 1:H
                xi(h, h_prime, h_bold) = ((h >= h_bold + 1) - sum(s_h(h_bold + 1:H))) * ((h_bold >= h_prime) - sum(s_h(1:h_bold)));
                temp(h, h_prime, h_bold) = xi(h, h_prime, h_bold) * (1 / rho_h(h_bold));
            end
            if h < h_prime
                epsilon_h_compl(h, h_prime) = sum(temp(h, h_prime, 1:H-1));
            end
        end
    end
end
