
function mpl = margProdLabor2(theta, kappa, z, alphaVec, xT)

    % note: This is a less efficient version. If you already have calculated the labor demand corresponding to
    % these inputs, it is more efficient to use the alternative version of this
    % function that takes labor demands as inputs.
    % Inputs:
    %   theta - blueprint scale (scalar)
    %   kappa - blueprint shape (scalar)
    %   z - productivity (scalar)
    %   alphaVec - array of comparative advantage values (vector)
    %   xT - array of H-1 task thresholds (vector)
    %
    % Output:
    %   mpl - marginal products of labor (vector)
    
    l = TaskBasedProduction.unitInputDemand(theta, kappa, z, alphaVec, xT);
    mpl = TaskBasedProduction.margProdLabor(l, alphaVec, xT);
end