clc;
clear all;
close all;

data = readtable('conductivity_data.csv'); % Assuming the table has columns 'phi' and 'sigma'

% conversion wt to volume fraction
% assuming values of density for CB and UHMWPE
ro_CB = 1.9; %g/cm^3
ro_UHMWPE = 0.93; %g/cm^3
% get rid of percents
data.phi = data.phi / 100;
phi = data.phi ./ (data.phi + (1 - data.phi) .* (ro_CB ./ ro_UHMWPE));
% phi = data.phi / 100.0;
sigma = data.sigma;

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
plot(phi, sigma, '^','MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'DisplayName', 'Data');
hold on;

% Plot the implicit curve
fimplicit(implicit_expr_sub, [0 0.2 0 1e3]); % Adjust the range as needed
set(gca, 'YScale', 'log'); % Set y-axis to logarithmic scale
%title(['Data and Fitted Implicit Curve, R^2 = ', num2str(r_squared)]);
xlabel('V_{CB}');
ylabel('\sigma [S/m]');
legend('Experiments', 'Effective medium fit', 'Location', 'best');
ax = gca;
ax.FontSize = 12;


% Display the fitted parameters
disp('Fitted parameters:');
disp(['t = ', num2str(params(1))]);
disp(['phi_c = ', num2str(params(2))]);

%% finding values for a fitted function
% values for weight fractions 0.5, 1, 1.5, 2.5, 5, 10, 15%

% layer values
% phi_values = [0.022;0.044;0.067;0.112;0.226;0.465;0.716];

% composite values
% phi_values = [0; 0.0025; 0.0049; 0.0074; 0.0124; 0.0251; 0.0516; 0.0795];
phi_values = phi;

L = length(phi_values);
sigma_values = zeros(L,1);

for i=1:L
    substituted_expr = subs(implicit_expr_sub, phi_sym, phi_values(i));
    sigma_values(i) = double(solve(substituted_expr, sigma_sym));
end

model_output = table(phi_values, sigma_values);
model_output.Properties.VariableNames = {'phi', 'sigma'};
model_output
%%
figure(3)
plot(phi, sigma, '^', 'MarkerSize', 7, 'MarkerFaceColor', '#33a2ff', 'MarkerEdgeColor', '#33a2ff')
set(gca, 'YScale', 'log');
% title("Experimental data only")
xlabel('V_{CB}');
ylabel('\sigma [S/m]');
ylim([0 1e3])
ax = gca;
ax.FontSize = 14;

%% extracting values from FEA results
data_FEA = readtable('auto_properties_thermal.xlsx');
data_FEA.wt__1 = data_FEA.wt__1 / 100; % no percents
phi_FEA = data_FEA.wt__1 ./ (data_FEA.wt__1 + (1 - data_FEA.wt__1) .* (ro_CB ./ ro_UHMWPE));
%add zero at the beg. 
phi_FEA = [0;phi_FEA];
sigma_FEA = [1e-10; data_FEA.K_avg_1];
%% experimental / medium / FEA results comparison
% experimental dataset (phi,sigma)
% curve fit dataset (phi_values, sigma_values)
% FEA dataset (phi_FEA, K_avg_1 (column name from file))

figure(4)
plot(phi([1:end-1]), sigma([1:end-1]), 'v', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b')
hold on
plot(phi_values([1:end-1]), sigma_values([1:end-1]), '*', 'MarkerFaceColor', '#FFA500', 'MarkerEdgeColor', '#FFA500')
plot(phi_FEA([1:end-1]), sigma_FEA([1:end-1]), 'o-', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'Color', 'g')
hold off
set(gca, 'YScale', 'log');
xticks([0:0.01:max(phi), phi(end)])
ylim([0 10^3])
xlabel('V_{CB}');
ylabel('\sigma [S/m]');
legend("Experimental", "Effective medium fit", "FEA results", 'Location',...
    'best', 'Box', 'off');
ax = gca;
ax.FontSize = 12;

%% check if K_eff is linear with K_layer

% values from effective medium equation
K_layer = [1e-10; 0.80; 4.71; 13.35; 46.65; 252.24; 1409.96; 3934.13];

% check linearity with linear fit curve on data
p = polyfit(K_layer, sigma_FEA, 1); % Linear fit (degree 1)
yfit = polyval(p, K_layer); % Evaluate the fit
yresid = sigma_FEA - yfit; % Residuals
SSresid = sum(yresid.^2); % Residual sum of squares
SStotal = (length(sigma_FEA)-1) * var(sigma_FEA); % Total sum of squares
rsq = 1 - SSresid/SStotal; % R² calculation

% coefficients for Kl vs. K_eff
p_r = polyfit(sigma_FEA, K_layer, 1);
% Display the equation for reference
fprintf('The linear equation is: K_l = %.2f*K_FEA + %.6f\n', p_r(1), ...
    p_r(2));

% plot
figure(5)
plot(K_layer, sigma_FEA, 'o') % original
hold on
plot(K_layer, yfit, '-') % Line of fit
xlabel('K_l')
ylabel('\sigma_F')
title(['Linear Fit with R2 = ', num2str(rsq),...
    '. Slope is ', num2str(round(p(1), 2)), '; Intercept is ', num2str(p(2))]);
legend("Experimental data", "Linear fit", 'Location', 'northwest')
hold off

%% now fitting FEA curve to experimental

% first, let's extract FEA points that coincide on vf with experimental
% I just get positions visually from the previous graph and reassign to
% another variable here

% experimental
phi_fit = phi([2,4,5,6,8]);
sigma_fit = sigma([2,4,5,6,8]);

% FEA
phi_FEA_fit = phi_FEA([2,3,5,6,7]);
sigma_FEA_fit = sigma_FEA([2,3,5,6,7]);
% (should be phi_fit == phi_FEA_fit)

% now to find what K_l should be in reality, using K_l vs K_Eff coef.
K_layer_eval = polyval(p_r, sigma_fit)

%% find a "real" curve based on adjusted K_layer values
% K adjusted for reverse fit
phi_adj = [0.022, 0.044, 0.112, 0.226, 0.465]'; % vf carbon in layer
phi_adj_wt = wtp_to_vf(phi_adj, 1);
sigma_adj = [0.8; 4.71; 46.65; 252.24; 1409.96];

% reverse fit
[expression_adj, params_adj, phi_adj] = eff_medium(phi_adj_wt, sigma_adj);

%% plot initial curve vs. adjusted
figure(6)
fimplicit(implicit_expr_sub, [0 0.2 1e-10 1e4], 'LineWidth', 1.5, ...
    'Color', '#ff8d33'); 
hold on
fimplicit(expression_adj, [0 0.2 1e-10 1e4], 'LineWidth', 1.5, ...
    'Color', '#33a2ff'); 
%plot(phi, sigma, 'v', 'MarkerFaceColor', '#ff8d33', 'MarkerEdgeColor', '#ff8d33')
%plot(phi_adj, sigma_adj, 'v', 'MarkerFaceColor', '#33a2ff', 'MarkerEdgeColor', '#33a2ff')
hold off
set(gca, 'YScale', 'log'); % Set y-axis to logarithmic scale
%title(['Data and Fitted Implicit Curve, R^2 = ', num2str(r_squared)]);
xlabel('V_{CB}');
ylabel('\sigma [S/m]');
legend(['Random CB distribution | \phi_c = ', num2str(round(params(2),4)),...
    ' | t = ', num2str(round(params(1),2))],...
    ['\muCT-driven CB distribution | \phi_c = ', num2str(round(params_adj(2),4)),...
    ' | t = ', num2str(round(params_adj(1),2))], 'location', 'best');
ax = gca;
ax.FontSize = 12;

figure(7)
fimplicit(implicit_expr_sub, [0 0.5 1e-10 1e4], 'LineWidth', 1.5, ...
    'Color', '#ff8d33'); 
hold on
fimplicit(expression_adj, [0 0.5 1e-10 1e4], 'LineWidth', 1.5, ...
    'Color', '#33a2ff'); 
plot(phi, sigma, 'v', 'MarkerFaceColor', '#ff8d33', 'MarkerEdgeColor', '#ff8d33')
plot(phi_adj, sigma_adj, 'v', 'MarkerFaceColor', '#33a2ff', 'MarkerEdgeColor', '#33a2ff')
hold off
set(gca, 'YScale', 'log'); % Set y-axis to logarithmic scale
%title(['Data and Fitted Implicit Curve, R^2 = ', num2str(r_squared)]);
xlabel('V_{CB}');
ylabel('\sigma [S/m]');
legend(['Random CB distribution | \phi_c = ', num2str(round(params(2),4)),...
    ' | t = ', num2str(round(params(1),2))],...
    ['\muCT-driven CB distribution | \phi_c = ', num2str(round(params_adj(2),4)),...
    ' | t = ', num2str(round(params_adj(1),2))], 'location', 'best');
ax = gca;
ax.FontSize = 12;







