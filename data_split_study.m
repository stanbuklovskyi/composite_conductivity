clc; clear;

%% import experimental data
data = readtable('conductivity_data.csv'); % Assuming the table has columns 'phi' and 'sigma'

%% set up random split 
phi = data.phi;
sigma = data.sigma;

% how many times to do the stduy
numSubsamples = 1;
for i = 1:numSubsamples
    
    % number of points to leave
    numPoints = randi([2,length(phi) - 4]);
    
    indices = randperm(length(phi) - 1, numPoints);
    phi_train = sort(phi(indices))
    sigma_train = sort(sigma(indices))
    
    % fit the curve using custom function
    [expression, params, phi] = eff_medium(phi_train,sigma_train);
end

%% saving data for future plots
% data_out = cell(6,5);
data_out{6,1} = phi_train; 
data_out{6,2} = sigma_train;
data_out{6,3} = expression;
data_out{6,4} = params(1); % t parameter
data_out{6,5} = params(2); % critical phi parameter