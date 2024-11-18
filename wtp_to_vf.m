function [out_array] = wtp_to_vf(in_array, reverse, ro_p, ro_m)
% wtp_to_vf Function Summary
% This function transforms wt% to volume fraction...
% or the other way around if reverse = True
%
% Input:
%
% "in_array": array of values that would be recalculated...
% to vf or wt% depending on reverse input
%
% "reverse": flag that if set to true calculates...
% wt% from volume fraction
%
% "ro_p": density of particle material
% "ro_m": density of matrix material
% Note: ro_p and ro_m should be the same dimensions.
%
%----------------------------------------------------------

% set up default value of reverse as False
arguments
    in_array;
    reverse = 0;
    % default are for CB and UHMWPE
    ro_p = 1.9;
    ro_m = 0.93;
end

if ~reverse
    in_array = in_array / 100.0;
    out_array = in_array ./ (in_array + (1 - in_array) .* (ro_p ./ ro_m));
else
    disp("This part doesn't work now")
end


end

