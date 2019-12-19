%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get the input characteristic of a pattern
% 
% Created: 08/31/2018 (mm/dd/yyyy)
% Modfied: 02/26/2019 (mm/dd/yyyy) Hexagon pattern
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #y_c, x_c# the position of the pattern
%
% #amp# the amplitude of the pattern
%
% #offset# the offset of the pattern
%
% #flag_pattern# flag to specify the pattern to fit
%   'Gaussian'
%   'Hexagon'
%   'Moffat'
%
% #par_pat# parameters of the pattern
% 
%%%%%%%%
% Ouput
%%%%%%%%
% #sel_var# the selected variable
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pattern = get_pattern(y_c, x_c, amp, offset, flag_pattern, ...
    par_pat)
    switch flag_pattern
        case {'Gaussian', 'Hexagon'}
            % Input: y_c, x_c, sig, amp, offset
            if numel(par_pat)==0
                pattern = {y_c, x_c, [], amp, offset} ;
            else
                pattern = {y_c, x_c, par_pat(:,:,1), amp, offset} ;
            end
        case 'Moffat'
            % Input: y_c, x_c, alpha, beta, amp, offset
            if numel(par_pat)==0
                pattern = {y_c, x_c, [], [], amp, offset} ;
            else
                pattern = {y_c, x_c, par_pat(:,:,1), par_pat(:,:,2), ...
                    amp, offset} ;
            end
        otherwise
            error([flag_pattern, ...
                ': Unkown model for the pattern...']) ;
    end
end