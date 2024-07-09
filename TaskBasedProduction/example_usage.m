clear
clc
% Initialize parameters
theta = 1.0;
kappa = 0.5;
z = 1.2;
alphaVec = [0.1, 0.2, 0.3];
xT = [0.4, 0.5];

% Call unitInputDemand and print the output
labor_demand = TaskBasedProduction.unitInputDemand(theta, kappa, z, alphaVec, xT);
disp('Labor Demand:');
disp(labor_demand);

% Call margProdLabor with labor demand
mpl = TaskBasedProduction.margProdLabor(labor_demand, alphaVec, xT);
disp('Marginal Products of Labor:');
disp(mpl);

% Call margProdLabor with blueprint characteristics
mpl_theta = TaskBasedProduction.margProdLabor2(theta, kappa, z, alphaVec, xT);
disp('Marginal Products of Labor (using blueprint characteristics):');
disp(mpl_theta);

% Call prod_fun with labor demand
[q, xT] = TaskBasedProduction.prod_fun(labor_demand, theta, kappa, z, alphaVec);
disp('Quantity Produced:');
disp(q);
disp('Task Thresholds:');
disp(xT);

% Call elasticity_sub_comp
[epsilon_h_sub, epsilon_h_compl] = TaskBasedProduction.elasticity_sub_comp(xT, labor_demand, q, mpl, theta, kappa, z, alphaVec);

disp('Elasticity of Substitution:');
disp(epsilon_h_sub);
disp('Elasticity of Complementarity:');
disp(epsilon_h_compl);

%% General parameterization

% Define b_g as a function handle
b_g = @(x) (x.^(kappa-1) .* exp(-x/theta)) / (theta^kappa * gamma(kappa));

% Define e_h functions as function handles
e_h1 = @(x) exp(0.1 * x);
e_h2 = @(x) exp(0.2 * x);
e_h3 = @(x) exp(0.3 * x);
e_h = {e_h1, e_h2, e_h3}; % Cell array of function handles

% Call the unitInputDemand_general function
labor_demand_general = TaskBasedProduction.unitInputDemand_general(xT, z, b_g, e_h);

% Display the result
disp('Labor Demand:');
disp(labor_demand_general);

[q_gen, xT_gen, fval]= TaskBasedProduction.prod_fun_general(labor_demand_general,z,b_g, e_h);
mpl_gen=TaskBasedProduction.margProdLabor_general(xT_gen, labor_demand_general, e_h);

[epsilon_h_sub_gen, epsilon_h_compl_gen] = TaskBasedProduction.elasticity_sub_comp_general(xT_gen, labor_demand_general, q_gen, mpl_gen, z, b_g, e_h);

