clc; clear;

%% test

% examplary file
filename = 'example_model_0';
new_property_number = 264;

new_file = change_material_properties(filename,new_property_number);

%% Linux bash commands

% 3 thermal loadcases
loadcases = ['example_model_0', 'example_model_0', 'example_model_0'];

new_property_number = 364;

for i = 1:length(loadcases)
    filename = loadcases(i);
    new_file = change_material_properties(filename,new_property_number);
    marc_run = ['marc2019 -jid ', new_file];
    [status, ~] = system(marc_run);

    if status == 0
        disp(['File ', new_file, ' was ran in Marc succesfully.'])
    else
        disp(['Error! File ', new_file, ' was not ran in Marc succesfully.'])
    end
    
end

%% testing
Marc_format(240)
