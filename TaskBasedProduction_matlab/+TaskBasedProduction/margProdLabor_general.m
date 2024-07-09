function mpl = margProdLabor_general(xT, l, e_h)
    % Calculate the marginal products of labor for each worker type
    % Inputs:
    %   xT - array of H-1 task thresholds
    %   l - array of H labor demand values
    %   e_h - cell array of comparative advantage functions e_h{h}

    H = length(e_h);
    mpl_over_mpl1 = 1.0;  % Initialize with 1.0 for h=1
    temp = zeros(H-1, 1); % Pre-allocate array for ratio values

    % Calculate the ratio e_{h} / e_{h-1} for h = 2:H and evaluate at xT[h-1]
    for h = 2:H
        ratio_value = e_h{h}(xT(h-1)) / e_h{h-1}(xT(h-1));
        temp(h-1) = ratio_value;
    end

    % Calculate cumulative product for mpl_over_mpl1
    mpl_over_mpl1 = [1; cumprod(temp)];

    % Calculate mpl1
    mpl1 = 1 / sum(mpl_over_mpl1 .* l);

    % Calculate mpl
    mpl = mpl_over_mpl1 * mpl1;
end
