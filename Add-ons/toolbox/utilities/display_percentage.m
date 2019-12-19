%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to display the evolution of the percentage of a loop in the
% command window
%
% Created: 09/03/2018 (mm/dd/yyyy)
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #flag_step# flag to specify the step of the percentage
%   'init' -> Initialization / varin = text
%   'iter' -> Iteration / varin = {per, per_disp, delta_per, lastprint}
%   'exit' -> Ending iteration / varin = lastprint
%
% #text# text to describe the percentage
%
% #per# current value of the percentage
%
% #per_disp# currently displayed percentage
%
% #delta_per# incremental value at which the percentage is displayed
%
% #lastprint# last printed message
%
%%%%%%%%
% Ouput
%%%%%%%%
% #per_disp_out# new displayed percentage
%
% #lastprint_out# new last printed message
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [per_disp_out, lastprint_out] = display_percentage(flag_step, ...
    varin)

    switch flag_step
        case 'init' % varin = text
            fprintf(varin) ;
            lastprint_out = fprintf('...') ;
            per_disp_out = 0 ;
        case 'iter' % {per, per_disp, delta_per, lastprint}
            if varin{1}>varin{2}
                fprintf(repmat('\b', 1, varin{4}));
                lastprint_out = fprintf([': ', num2str(varin{2}), ' %%']) ;
                per_disp_out = (floor(varin{1}/varin{3})+1)*varin{3} ;
            else
                lastprint_out = varin{4} ;
                per_disp_out = varin{2} ;
            end
        case 'exit' % varin = lastprint
            fprintf(repmat('\b', 1, varin)) ;
            per_disp_out = fprintf(': done!') ;
            lastprint_out = [] ;
            disp(' ') ;
        otherwise
            error([flag_step, ...
                ': wrong possibility for ''flag_step''...']) ;
    end
end