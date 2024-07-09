function mpl = margProdLabor(inputDemand, alphaVec, xT)
    % margProdLabor Calculates marginal products of labor
    %
    % Calculate marginal products of labor for each worker type given an array
    % of H labor demand values, an array of comparative advantage values alpha,
    % and an array of H-1 task thresholds xT that corresponds to that labor demand.
    %
    % Inputs:
    %   inputDemand - array of H labor demand values (vector)
    %   alphaVec - array of comparative advantage values (vector)
    %   xT - array of H-1 task thresholds (vector)
    %
    % Output:
    %   mpl - marginal products of labor (vector)
    
    mpl_over_mpl1 = [1 cumprod(exp(diff(alphaVec) .* xT))]';
    mpl1 = 1 / sum(mpl_over_mpl1 .* inputDemand);
    mpl = mpl_over_mpl1 * mpl1;
end

