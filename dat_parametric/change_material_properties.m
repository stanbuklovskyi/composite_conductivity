function [new_file] = change_material_properties(filename,new_property_number, search_word)

% change_material_properties: Function Summary
% This function takes some example .dat file "filename", changes material
% property number and saves new file as "new_file"
%
% Inputs:
%
% filename - string, name of the example file without '.dat' extension
% new_property_number - material property number you would like to have
% search_word - char, key word for the name of the material to look for.
% default "search_word" is 'interphase'
%
% Outputs: 
% new_file - name of the newly created .dat file with a property changed...
% to new_property_value. File is created in the current directory.
%-----------------------------------------------------------------------------------------------

% to do: 
% 1) Change setup for actual interphase key word
% 2) Change textscan read format to read some positions char (see .dat)
% 3) Don't forget to use Marc_format(new_property_number) for Marc format
%
%------------------------------------------------------------------------

arguments
    filename;
    new_property_number;
    search_word = 'interphase';
end

% read/write files
file_id = fopen([filename, '.dat'], 'r');

new_file = [filename, '_prop', num2str(new_property_number), '.dat'];
new_file_id = fopen(new_file, 'w');

while ~feof(file_id)

    line = fgetl(file_id);

    if contains(line, search_word)
        % write it to a new file too
        fprintf(new_file_id, '%s\n', line);

        % skip a line with writting it
        line = fgetl(file_id);
        fprintf(new_file_id, '%s\n', line);
        
        % read and change property
        line_with_numbers = fgetl(file_id);
        numbers = textscan(line_with_numbers, '%10f%10f%10f%10f%10f%10f');
        numbers{3} = new_property_number;

        % Create the modified line
        modified_line = sprintf('%10.0f%10.0f%10.0f%10.0f%10.0f%10.0f', numbers{1},...
            numbers{2}, numbers{3}, numbers{4}, numbers{5}, numbers{6});
        
        % Write the modified line to the new file
        fprintf(new_file_id, '%s\n', modified_line);
    else
        % Write the line to a new file
        fprintf(new_file_id, '%s\n', line);
    end

end

% close files
fclose(file_id);
fclose(new_file_id);

end