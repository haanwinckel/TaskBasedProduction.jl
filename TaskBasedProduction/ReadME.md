# TaskBasedProduction

TaskBasedProduction is a MATLAB toolbox that provides functions for calculating unit labor demands, marginal products of labor, production function, assignment thresholds, elasticities of substitution and complementarities among worker types in a task-based production model. The package includes utilities for handling incomplete gamma functions and power series representations to facilitate these calculations.

## Installation

To install TaskBasedProduction, you can clone the repository and run the `install.m` script to add it to your MATLAB path:
## Usage Example
```matlab

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

% Call prod_fun and print the output
[q, xbar] = TaskBasedProduction.prod_fun(labor_demand, theta, kappa, z, alphaVec);
disp('Quantity Produced:');
disp(q);
disp('Task Thresholds:');
disp(xbar);

% Call margProdLabor with labor demand and print the output
mpl = TaskBasedProduction.margProdLabor(labor_demand, alphaVec, xT);
disp('Marginal Products of Labor (with labor demand):');
disp(mpl);

% Call margProdLabor with blueprint characteristics and print the output
mpl_alt = TaskBasedProduction.margProdLabor_theta(theta, kappa, z, alphaVec, xT);
disp('Marginal Products of Labor (with blueprint characteristics):');
disp(mpl_alt);

% Call elasticity_sub_comp with labor demand, MPL, xT and parameters of the gamma function and print the two outputs
epsilon_sub, epsilon_compl = TaskBasedProduction.elasticity_sub_comp(xT, labor_demand, q, mpl, theta, kappa, z, alphaVec);
disp('Allen partial elasticity of substitution:');
disp(epsilon_sub);
disp('Hicks partial elasticity of substitution:');
disp(epsilon_compl);
```
## Functions and Features
1) unitInputDemand: Calculates unit labor demands given blueprint scale theta, blueprint shape kappa, productivity z, an array of comparative advantage values alphaVec, and an array xT of thresholds in task space.
unitInputDemand(theta, kappa, z, alphaVec, xT, skipParamChecks)

2) margProdLabor: Calculates marginal products of labor for each worker type given an array of labor demands (margProdLabor) or given blueprint characteristics (margProdLabor2). Note that if the labor demand has been already calculated, the first function is more efficient

margProdLabor(inputDemand, alphaVec, xT)
margProdLabor2(theta, kappa, z, alphaVec, xT)



3) prod_fun: Calculates the quantity produced (q) and task thresholds (xT) given labor inputs (l), blueprint scale theta, blueprint shape kappa, productivity z, and an array of comparative advantage values alphaVec with H elements (one for each worker type).

prod_fun(l, theta, kappa, z, alphaVec)

4) elasticity_sub_comp: Calculates the elasticity of substitution and complementarities for a given labor inputs (l), an array xT of thresholds in task space (dimension H-1), array of marginal product of labor (MPL) for each H labor types, blueprint scale theta, blueprint shape kappa, productivity z, and an array of comparative advantage values alphaVec with H elements (one for each worker type). The function returns two matrices representing the elasticity of substitution and complementarity values for each worker type h (rows) relative to worker type h_prime (columns).

elasticity_sub_comp(xT, l, q, MPL, theta, kappa, z, alphaVec)

5) unitInputDemand_general
Calculates unit labor demands given an array xT of H-1 thresholds in task space, productivity value z, 
a density function b_g for the task distribution, and an array e_h of H functions
representing the cost of each labor type as a function of task complexity.

unitInputDemand_general(xT, z, b_g, e_h);

6)  margProdLabor_general
Calculate marginal products of labor for each worker type given an array
of H labor demand values, a vector of comparative advantage functions e_h, and
an array of H-1 task thresholds xT that corresponds to that labor demand.

margProdLabor_general(xT_gen, labor_demand_general, e_h);


7) prod_fun_general

Calculates the quantity produced (q), and task thresholds (xbar)
given labor inputs (l),  productivity z, and general blueprint density function and a vector of efficiency functions, one for each labor type.

prod_fun_general(labor_demand_general,z,b_g, e_h);

8) elasticity_sub_comp_general
Calculates the elasticity of substitution and complementarities for a given labor inputs (l), an array xT of thresholds in task space (dimension H-1), array of marginal product of labor (MPL) for each H labor types, blueprint density function b_g, firm productivity z, and a vector of comparative advantage function e_h. The function returns two matrices representing the elasticity of substitution and complementarity values for each worker type h (rows) relative to worker type h_prime (columns).

elasticity_sub_comp_general(xT_gen, labor_demand_general, q_gen, mpl_gen, z, b_g, e_h);


## Contributing
Contributions are welcome! Please feel free to submit a pull request or open an issue if you have any suggestions or find any bugs.
