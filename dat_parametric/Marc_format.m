function formattedStr = Marc_format(num)
% Marc_format Function Summary
%
% This function takes any number and converts it in a specific format...
% for Marc .dat files. Example: 240 -> 2.40e+2 -> ' 2.400000000000000+2'
%
%
% Input:
% num: number, random regular matlab number to be transformed.
% 
% Output:
% formattedStr: string, converted input number in a marc format.
%
%
%-----------------------------------------------------------------------

    % Step 1: Convert the number to scientific notation with 15 decimal places
    sciStr = sprintf('%.15e', num);
    
    % Step 2: Replace 'e' with '+' or '-' and ensure the exponent is formatted correctly
    sciStr = strrep(sciStr, 'e+', '+');
    sciStr = strrep(sciStr, 'e-', '-');
    
    % Step 3: Ensure the exponent is in the correct format
    % Remove any leading zeros in the exponent
    sciStr = regexprep(sciStr, '([+-])0+(\d)', '$1$2');
    
    % Step 4: Adjust the format to match the desired output
    % Ensure the total length is 19 characters for the number and 1 space before it
    formattedStr = [' ' sciStr];
    
    % If the total length is less than 20, pad with spaces
    if length(formattedStr) < 20
        formattedStr = [formattedStr, repmat('0', 1, 20 - length(formattedStr))];
    end
end