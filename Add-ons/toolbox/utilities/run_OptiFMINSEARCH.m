%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to run a fminsearch optimization
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
% solution (must be empty for this method which is unconstrained)
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
function [par_out] = run_OptiFMINSEARCH(C, const_min, const_max, ...
    option_opti, par_in)


    %% Initialization
    if ~isempty(const_min) || ~isempty(const_max)
        warning(['Constraints are given whereas ''fminsearch'' is an ', ...
            'unconstrained solver... Constraints not taken into', ...
            ' account...']) ;
    end

    %% Optimizer declaration
    options = optimset() ;
    if isfield(option_opti, 'verbose')
        if option_opti.verbose
            options.Display = 'iter' ;
        else
            options.Display = 'none' ; 
        end
    else
        options.Display = 'none' ; 
    end
    if isfield(option_opti, 'maxiter')
        options.MaxIter = option_opti.maxiter ; 
    end 
    if isfield(option_opti, 'TolFun')
        options.TolFun = option_opti.TolFun ; 
    end 
    if isfield(option_opti, 'TolX')
        options.TolX = option_opti.TolX ; 
    end
    if isfield(option_opti, 'MaxFunEvals')
        options.MaxFunEvals = option_opti.MaxFunEvals ; 
    end
    
    %% Optimization
    if isa(C,'Cost')
        % GlobalBioIm formalism
        if isfield(option_opti, 'memoizeOpts')
            if option_opti.memoizeOpts
                C.memoizeOpts.apply = true ;
            end
        end
        par_out = fminsearch(@(x)(calculateCost(C, x)), par_in, options) ;
    else
        par_out = fminsearch(C, par_in, options) ;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to define the cost function and its gradient with the
% GlobalBioIm formalism
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #C# the cost function
%
% #x# the point at each the function must be computed
%
%%%%%%%%
% Ouput
%%%%%%%%
% #cost# value of the cost
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cost] = calculateCost(C, x)
    cost = C*x ;
end