function [epsilon_h_sub, epsilon_h_compl] = elasticity_sub_comp(xT, l, q, MPL, theta, kappa, z, alphaVec)
    % elasticity_sub_comp Calculates the elasticity of substitution and complementarity
    %
    % Calculates the elasticity of substitution and complementarity for a given set of parameters.
    %
    % Inputs:
    %   xT - Array of task thresholds with H-1 elements (vector)
    %   l - Array of labor inputs of different types with H elements (vector)
    %   q - Total output (scalar)
    %   MPL - Array of marginal products of labor for each worker type with H elements (vector)
    %   theta - Blueprint scale parameter (scalar or array indexed by h with H elements)
    %   kappa - Blueprint shape parameter (scalar or array indexed by h with H elements)
    %   z - Productivity parameter (scalar or array indexed by h with H elements)
    %   alphaVec - Array of comparative advantage values with H elements (vector)
    %
    % Returns:
    %   epsilon_h_sub - Matrix of elasticity of substitution values for each worker type h (rows) relative to worker type h_prime (columns)
    %   epsilon_h_compl - Matrix of elasticity of complementarity values for each worker type h (rows) relative to worker type h_prime (columns)

    H = length(alphaVec);
    rho_h = zeros(H, 1);
    s_h = (MPL .* l) / q;
    epsilon_h_sub = zeros(H, H);
    epsilon_h_compl = zeros(H, H);
    xT = [xT(:); Inf];  % Add highest thresholds for highest worker type
    e_h_T = exp(alphaVec .* xT);  % Denominator of rho_h

    for h = 1:H-1
        b_g_T = xT(h)^(kappa - 1) * (1 / (z * gamma(kappa) * theta^kappa)) * exp(-xT(h) / theta);
        rho_h(h) = b_g_T * MPL(h) * (1 / e_h_T(h)) * (1 / (alphaVec(h+1) - alphaVec(h)));
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
                xi(h, h_prime, h_bold) = ((h >= h_bold + 1) - sum(s_h(h_bold+1:end))) * ((h_bold >= h_prime) - sum(s_h(1:h_bold)));
                temp(h, h_prime, h_bold) = xi(h, h_prime, h_bold) * (1 / rho_h(h_bold));
            end
             if h < h_prime
                epsilon_h_compl(h, h_prime) = sum(temp(h, h_prime, 1:H-1));
            end
        end
    end
end
