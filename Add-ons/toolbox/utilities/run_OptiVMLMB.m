%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to run a VMLMB optimization
% 
% Created: 08/03/2018 (mm/dd/yyyy)
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
%
% #par_in# initial parameters
%
%%%%%%%%
% Ouput
%%%%%%%%
% #par_out# output parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [par_out] = run_OptiVMLMB(C, const_min, const_max, ...
    option_opti, par_in)

    %% Initialization
    if isfield(option_opti, 'memoizeOpts')
        if option_opti.memoizeOpts
            C.memoizeOpts.apply = true ;
        end
    end

    %% Optimizer declaration
    VMLMB=OptiVMLMB(C, const_min, const_max) ; 
    if isfield(option_opti, 'ItUpOut')
        VMLMB.ItUpOut = option_opti.ItUpOut ; 
    else
        VMLMB.ItUpOut = 1 ; 
    end 
    if isfield(option_opti, 'verbose')
        VMLMB.verbose = option_opti.verbose ; 
    else
        VMLMB.verbose = false ; 
    end 
    if isfield(option_opti, 'm')
        VMLMB.m = option_opti.m ; 
    else
        VMLMB.m = 3 ; 
    end 
    if isfield(option_opti, 'maxiter')
        VMLMB.maxiter = option_opti.maxiter ; 
    end 

    %% Optimization
    VMLMB.run(par_in);
    if ~isempty(VMLMB.xopt)
        par_in = VMLMB.xopt ;
    end
    par_out = par_in ;
end