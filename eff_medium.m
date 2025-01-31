function [implicit_expr_sub, params, phi] = eff_medium(data_phi,data_sigma)

%"eff_medium" function summary 
% This function is used to fit effective medium equation to experimental...
% ...conductivity data. See (Blaszkiewicz, McLachlan 1992, Polymer Eng Sci)
%
% Inputs:
% (data_phi, data_sigma) - 1D arrays of % wt fractions and
% conductivivities (S/m) respectively. There is wt% to volume fraction
% conversion below.
%
% Outputs: 
% 'params' - two element 1D array of fitted coefficients.
%
% params(1) - t exponent. params(2) - phi_c critical volume fraction... 
% ... for percolation.
%
% 'implicit_expr_sub' - symbolic expression of the equation from the
% referenced paper with t and phi_c found and substituted for plotting.
% Note: if symbolic expression is used, syms phi_sym and sigma_sym should
% be defined with "syms phi_sym sigma_sym" expression in the main code.
% Else, symbols are unrecognized
%
% 'phi' - volume fraction used and converted to vf from wt%.  
%--------------------------------------------------------------------------

% conversion wt% to volume fractions (using densities of constituents)
ro_CB = 1.9; %g/cm^3
ro_UHMWPE = 0.93; %g/cm^3
% get rid of percents
data_phi = data_phi / 100;
phi = data_phi ./ (data_phi + (1 - data_phi) .* (ro_CB ./ ro_UHMWPE));
sigma = data_sigma;

% add zeros at the beginning for numerical fit
%phi = [0;phi];
%sigma = [1e-10;sigma];

% Define constants
sigma_l = 1e-10;
sigma_h = 1 / (1.15e-4);

% Define the implicit function
implicit_func = @(params, phi, sigma) ...
    (1-phi).*(sigma_l.^(1/params(1)) - sigma.^(1/params(1))) ./ (sigma_l.^(1/params(1)) + (1-params(2)).*sigma.^(1/params(1))./params(2)) + ...
    phi.*(sigma_h.^(1/params(1)) - sigma.^(1/params(1))) ./ (sigma_h.^(1/params(1)) + (1-params(2)).*sigma.^(1/params(1))./params(2));

% Define the objective function for optimization (sum of squared residuals)
objective_func = @(params) sum((implicit_func(params, phi, sigma)).^2);

% Initial guess for parameters [t, phi_c]
initial_guess = [1, 0.002];

% Set bounds for the parameters
lower_bounds = [1, 0.002];
upper_bounds = [5, 0.5];

% Fit the model to the data using fminsearch
% opts = optimset('Display', 'off');
% params = fminsearch(objective_func, initial_guess, opts);
opts = optimoptions('fmincon', 'Display', 'off');
params = fmincon(objective_func, initial_guess, [], [], [], [], lower_bounds, upper_bounds, [], opts);

% Calculate residuals
residuals = implicit_func(params, phi, sigma);
ss_res = sum(residuals.^2);
ss_tot = sum((sigma - mean(sigma)).^2);
r_squared = 1 - (ss_res / ss_tot);

% Calculate residuals
% residuals = log(implicit_func(params, phi, sigma) + 1);
% ss_res = sum(residuals.^2);
% ss_tot = sum((log(sigma + 1) - mean(log(sigma + 1))).^2);
% r_squared = 1 - (ss_res / ss_tot);

% Define symbolic variables for plotting
syms phi_sym sigma_sym t phi_c

% Define the implicit expression for plotting
implicit_expr1 = (1-phi_sym)*(sigma_l^(1/t) - sigma_sym^(1/t)) / (sigma_l^(1/t) + (1-phi_c)*sigma_sym^(1/t) / phi_c);
implicit_expr2 = phi_sym*(sigma_h^(1/t) - sigma_sym^(1/t)) / (sigma_h^(1/t) + (1-phi_c)*sigma_sym^(1/t) / phi_c);
implicit_expr = implicit_expr1 + implicit_expr2;

% Substitute the fitted parameters into the implicit expression
implicit_expr_sub = subs(implicit_expr, {t, phi_c}, {params(1), params(2)});

% Plot the data points
figure;
plot(phi, sigma, 'v', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'DisplayName', 'Data');
hold on;

% Plot the implicit curve
fimplicit(implicit_expr_sub, [0 0.7 0 1e4]); % Adjust the range as needed
set(gca, 'YScale', 'log'); % Set y-axis to logarithmic scale
%title(['Data and Fitted Implicit Curve, R^2 = ', num2str(r_squared)]);
xlabel('V_{CB}');
ylabel('\sigma [S/m]');
legend('Experiments', ['Effective medium | \phi_c = ', num2str(round(params(2),4)), ...
    ' | t = ', num2str(round(params(1),2))], 'Location', 'best');
ax = gca;
ax.FontSize = 12;


% Display the fitted parameters
disp('Fitted parameters:');
disp(['t = ', num2str(params(1))]);
disp(['phi_c = ', num2str(params(2))]);
end

