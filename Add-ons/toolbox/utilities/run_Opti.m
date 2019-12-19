%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to run an optimization procedure
% 
% Created: 08/06/2018 (mm/dd/yyyy)
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #C# the cost function
%
% #const_min, const_max# minimal and maximal bound acceptable for the
% solution
%
% #option_opti# structure containing the options for the VMLMB optimizer
%   default: VMLMB
%
% #par_in# initial parameters
%
%%%%%%%%
% Ouput
%%%%%%%%
% #par_out# output parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [par_out] = run_Opti(C, const_min, const_max, ...
    option_opti, par_in)

    if ~isfield(option_opti, 'method')
       option_opti.method = 'VMLMB' ; 
    end

    switch option_opti.method
        case 'VMLMB'
            [par_out] = run_OptiVMLMB(C, const_min, const_max, ...
                option_opti, par_in) ;
        case 'fminunc'
            [par_out] = run_OptiFMINUNC(C, const_min, const_max, ...
                option_opti, par_in) ;
        case 'fminsearch'
            [par_out] = run_OptiFMINSEARCH(C, const_min, const_max, ...
                option_opti, par_in) ;
        otherwise
            error([option_opti.method, ...
                ':unknown optimization procedure...']) ;
    end
end