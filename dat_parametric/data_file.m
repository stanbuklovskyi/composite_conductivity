clc; clear;

%% setup of example data files

% 3 thermal loadcases
loadcases = ["composite_cond_T_xx_0_5", "composite_cond_T_yy_0_5", "composite_cond_T_zz_0_5"];
new_property_number = 364;

% check if files exist in the current folder
file_counter = 0;
for i = 1:length(loadcases)
    if isfile(strcat(loadcases(i),".dat"))
        file_counter = file_counter + 1;
        
        if file_counter == length(loadcases)
            disp('All .dat files exist')
        end
    else
       disp(strcat(strcat(loadcases(i),".dat"), " does not exist")); 
    end
end

%%

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
