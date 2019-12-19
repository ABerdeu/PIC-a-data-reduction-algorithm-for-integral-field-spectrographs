%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get the scaling factor of the residues for the reweighted
% least square optimization
% 
% Created: 07/26/2018 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #res# the current residues
%
% #flag_s# method to scale the residue before applying the
%   objective function (default: none)
%   ->  'none'	-> no scaling applied
%   ->  'MAD'   -> median absolute deviation
%   ->  'STD'   -> standard deviation
%   ->  numeric -> scaling
%
%%%%%%%%
% Ouput
%%%%%%%%
% #s# the scaling factor
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = get_scaling_residues_factor(res, flag_s)
    if ~isnumeric(flag_s)
        switch flag_s
            case 'STD'
                % standard deviation
                s = 3 * std(res(:)) ;
            case 'MAD'
                % median absolute deviation
                s = 1.4826 * median( ...
                    abs(res(:)-median(res(:)))) ;
            case 'none'
                s = 1 ;
            otherwise
                error(['Unknown possibility to', ...
                    ' compute the scaling factor...']) ;
        end
    else
        s = flag_s ;
    end
end