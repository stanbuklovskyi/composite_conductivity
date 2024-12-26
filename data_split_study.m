clc; clear;

%% import experimental data
data = readtable('conductivity_data.csv'); % Assuming the table has columns 'phi' and 'sigma'

%% set up random split 
phi = data.phi;
data_phi = data.phi / 100; % get rid of percents
ro_CB = 1.9; %g/cm^3
ro_UHMWPE = 0.93; %g/cm^3
% get colume fraction out of wt%
data_phi = data_phi ./ (data_phi + (1 - data_phi) .* (ro_CB ./ ro_UHMWPE));
sigma = data.sigma;

%%
% how many times to do the stduy
numSubsamples = 1;
for i = 1:numSubsamples
    
    % number of points to leave
    numPoints = randi([2,length(phi) - 1]);
    
    indices = randperm(length(phi) - 1, numPoints);
    phi_train = sort(phi(indices))
    sigma_train = sort(sigma(indices))
    
    % fit the curve using custom function
    [expression, params, phi] = eff_medium(phi_train,sigma_train);
end

%% saving data for future plots
% attention! This part of the code running is a manual process
% data_out = cell(8,5); % cell array to manually store the data for 8 cases
% set up as many cases as you want

% data_out{4,1} = phi; 
% data_out{4,2} = sigma_train;
% data_out{4,3} = expression;
% data_out{4,4} = params(1); % t parameter
% data_out{4,5} = params(2); % critical phi parameter

%% plots

% one plot figure has:
% Full experimental data - not filled markers
% Experimental data used for fit (data_out{i,1} vs. data_out{i,2})
% plot of the curve expression data_out{i,3} with parameters...
% data_out{i,4} for t, and data_out{i,5} for phi_c

% i is the i-th figure
% all figure are exported separately as "plot_(i).jpg" and...
% added to a gif for a quick movie

% this is needed to make sure symbols in function output expression are defined
syms phi_sym sigma_sym

figure;
% experimental
plot(data_phi, data.sigma, 'v', 'MarkerEdgeColor', 'b','Linewidth', 1.1, 'MarkerSize', 7, 'DisplayName', 'Expresimental');
hold on;

% training set
plot(data_out{6,1}, data_out{6,2}, 'v','MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b','MarkerSize', 7, 'DisplayName', 'Training');

% Plot the implicit curve
fimplicit(data_out{6,3}, [0 0.15 0 1e3], 'LineWidth', 1.4); % Adjust the range as needed
set(gca, 'YScale', 'log'); % Set y-axis to logarithmic scale
%title(['Data and Fitted Implicit Curve, R^2 = ', num2str(r_squared)]);
xlabel('\phi volume fraction');
ylabel('\sigma conductivity');
legend('Experiments', 'Training', ['Effective medium | \phi_c = ', num2str(round(params(2),4)), ...
    ' | t = ', num2str(round(params(1),2))], 'Location', 'best');
ax = gca;
ax.FontSize = 12;
yticks(logspace(-10,3,14))
xticks(linspace(0,0.15, 7))
ylim([0,1e3])
hold off

%% plots in a loop 
% output in one animated git + all figs saved separately as .jpg
syms phi_sym sigma_sym

% separate plots
for i = 1:size(data_out, 1)
    figure;
    % experimental
    plot(data_phi, data.sigma, 'v', 'MarkerEdgeColor', 'b','Linewidth', 1.1, 'MarkerSize', 7, 'DisplayName', 'Expresimental');
    hold on;
    
    % training set
    plot(data_out{i,1}, data_out{i,2}, 'v','MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b','MarkerSize', 7, 'DisplayName', 'Training');
    
    % Plot the implicit curve
    fimplicit(data_out{i,3}, [0 0.15 0 1e3], 'LineWidth', 1.4); % Adjust the range as needed
    set(gca, 'YScale', 'log'); % Set y-axis to logarithmic scale
    %title(['Data and Fitted Implicit Curve, R^2 = ', num2str(r_squared)]);
    xlabel('\phi volume fraction');
    ylabel('\sigma conductivity');
    legend('Experiments', 'Training', ['Effective medium | \phi_c = ', num2str(round(data_out{i,5},4)), ...
        ' | t = ', num2str(round(data_out{i,4},2))], 'Location', 'best');
    ax = gca;
    ax.FontSize = 12;
    yticks(logspace(-10,3,14))
    xticks(linspace(0,0.15, 7))
    ylim([0,1e3])
    hold off
    
    saveas(gcf, ['plot_', num2str(i), '.jpg']);
    close;
end

% create gif with previously created images
for i = 1:size(data_out, 1)
    img = imread(['plot_', num2str(i), '.jpg']);
    [imind, cm] = rgb2ind(img, 256);
    if i == 1
        imwrite(imind, cm, 'movie.gif', 'gif', 'Loopcount', inf, 'DelayTime', 3);
    else
        imwrite(imind, cm, 'movie.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 3);
    end
end
