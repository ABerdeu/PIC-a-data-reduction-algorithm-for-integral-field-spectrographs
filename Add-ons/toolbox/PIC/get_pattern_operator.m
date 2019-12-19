%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get a pattern operator
% 
% Created: 05/15/2018 (mm/dd/yyyy)
% Modified: 05/17/2018 (mm/dd/yyyy) More possibilities for the gaussian
% simulation
% Modified: 09/03/2018 (mm/dd/yyyy) moving from Gaussian pattern to user
% defined pattern
% Modfied: 02/26/2019 (mm/dd/yyyy) Hexagon pattern
% Modfied: 03/11/2019 (mm/dd/yyyy) Getting the operator
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #[list_y, list_x]# list of the position of the pattern in the figure
%
% #y_pat, x_pat# position of the pattern
% 
% #par_pat# parameters of the pattern
%
% #amp_pat# amplitude of the pattern
%
% #off_pat# offset of the pattern
%
% #pattern_model# flags to specify the pattern to fit
%   pattern_model.flag_norm -> Normalized?
%   pattern_model.flag_profile -> Gaussian / Hexagon / Moffat
%   pattern_model.oversampling -> oversampling of the pattern (default: 1)
%
%%%%%%%%
% Ouput
%%%%%%%%
% #op_pattern# the pattern operator depending on the parameters that are
% not defined
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function op_pattern = get_pattern_operator(list_y, list_x, y_pat, ...
    x_pat, par_pat, amp_pat, off_pat, pattern_model)

    % Oversampling the pattern if needed    
    if isfield(pattern_model, 'oversampling') && ...
            pattern_model.oversampling>1
        [OpDS, list_y_OS, list_x_OS] = get_downsampling_operator( ...
            pattern_model.oversampling, list_y, list_x) ;
    else
        OpDS = 1 ;
        list_y_OS = list_y ;
        list_x_OS = list_x ;
    end
    
    % Getting the pattern operator
    switch pattern_model.flag_profile
        case 'Gaussian'
            % Input: y_c, x_c, sig, amp, offset
            Pattern = OpGauss(list_y_OS, list_x_OS, y_pat, x_pat, ...
                par_pat, amp_pat, off_pat, pattern_model.flag_norm) ;
        case 'Hexagon'
            % Input: y_c, x_c, sig, amp, offset
            Pattern = OpPSF_Hex(list_y_OS, list_x_OS, y_pat, x_pat, ...
                pattern_model.theta, par_pat, amp_pat, off_pat) ;
        case 'Moffat'
            % Input: y_c, x_c, alpha, beta, amp, offset
            Pattern = OpMoffat(list_y_OS, list_x_OS, y_pat, x_pat, ...
                par_pat(1), par_pat(2), amp_pat, off_pat, ...
                pattern_model.flag_norm) ;
        otherwise
            error([pattern_model.flag_profile, ...
                ': Unkown model for the pattern...']) ;
    end

    % Final operator
    op_pattern = OpDS*Pattern ;
end