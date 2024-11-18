clc;clear;

%% data import 
data = readtable('conductivity_data.csv'); % Assuming the table has columns 'phi' and 'sigma'

%%
% curve fit by function
[implicit_expr_sub, params, phi] = eff_medium(data.phi, data.sigma); 

%% finding values for a fitted function
% values for weight fractions 0.5, 1, 1.5, 2.5, 5, 10, 15%

% layer values
phi_values = [0.022;0.044;0.067;0.112;0.226;0.465;0.716];

% composite values
% phi_values = [0; 0.0025; 0.0049; 0.0074; 0.0124; 0.0251; 0.0516; 0.0795];
% phi_values = phi;

L = length(phi_values);
sigma_values = zeros(L,1);
syms phi_sym sigma_sym

for i=1:L
    substituted_expr = subs(implicit_expr_sub, phi_sym, phi_values(i));
    sigma_values(i) = double(solve(substituted_expr, sigma_sym));
end

model_output = table(phi_values, sigma_values);
model_output.Properties.VariableNames = {'phi', 'sigma'};
model_output

%% reverse fit the curve 
model_output_phi = wtp_to_vf(model_output.phi, 1);
model_output_sigma = model_output.sigma;
[implicit_expr_sub_r, params_r, phi_r] = eff_medium(model_output_phi, model_output_sigma); 



