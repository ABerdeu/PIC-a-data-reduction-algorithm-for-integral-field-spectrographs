%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to run a fminunc optimization
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
function [par_out] = run_OptiFMINUNC(C, const_min, const_max, ...
    option_opti, par_in)


    %% Initialization
    if ~isempty(const_min) || ~isempty(const_max)
        warning(['Constraints are given whereas ''fminunc'' is an ', ...
            'unconstrained solver... Constraints not taken into', ...
            ' account...']) ;
    end
    if isfield(option_opti, 'memoizeOpts')
        if option_opti.memoizeOpts
            C.memoizeOpts.apply = true ;
        end
    end

    %% Optimizer declaration
    options = optimoptions(@fminunc, ...
        'Algorithm','quasi-newton', ...
        'GradObj', 'on') ;
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
        options.MaxIterations = option_opti.maxiter ; 
    end 

    %% Optimization
    if isa(C,'Cost')
        % GlobalBioIm formalism
        par_out = fminunc(@(x)(calculateCost(C, x)), par_in, options) ;
    else
        par_out = fminunc(C, par_in, options) ;
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
% #gradient# value of the cost's gradient
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cost, gradient] = calculateCost(C, x)
    cost = C*x ;
    gradient = C.applyGrad(x) ;
end