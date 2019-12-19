%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to compute the profile of a symmetric pattern
% 
% Created: 09/04/2018 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #pos# the list of the positions on which the profile must be computed
%   (list of 1D, cell if 2D)
%
% #amp# the amplitude of the profile
%
% #offset# the offset of the profile
%
% #pattern_model# flags to specify the pattern to fit
%   pattern_model.flag_norm -> Normalized?
%   pattern_model.flag_profile -> Gaussian / Moffat
%
% #par_pat# parameters of the pattern
%
% #rate_OS# oversampling the profile?
%
%%%%%%%%
% Ouput
%%%%%%%%
% #profile# the pattern profile
%
% #amplitude# the amplitude of the profile
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [profile, amplitude] = get_pattern_profile(pos, amp, offset, ...
    pattern_model, par_pat, rate_OS)

    if nargin > 5 && rate_OS > 1
        if iscell(pos)
            [OpDS, list_y, list_x] = get_downsampling_operator( ...
                rate_OS, pos{2}, pos{1}) ;
            [list_x, list_y] = meshgrid(list_x, list_y) ;
        else
            [OpDS, list_y] = get_downsampling_operator( ...
                rate_OS, pos) ;
            list_x = 0 ;
        end
        r = (list_y.^2+list_x.^2).^0.5 ;
    else
        OpDS = 1 ;
        if iscell(pos)
            [list_x, list_y] = meshgrid(pos{2}, pos{1}) ;
        else
            list_y = pos ;
            list_x = 0 ;
        end
        r = (list_y.^2+list_x.^2).^0.5 ;
    end

    switch pattern_model.flag_profile
        case 'Gaussian'
            if pattern_model.flag_norm
                profile = amp/(2*pi*par_pat(1)^2) .* ...
                    exp(-1/2*r.^2/par_pat(1)^2) + offset ;
                amplitude = amp/(2*pi*par_pat(1)^2) ;
            else
                profile = amp.*exp(-1/2*r.^2/par_pat(1)^2) + offset ; 
                amplitude = amp ;
            end
        case 'Moffat'
            if pattern_model.flag_norm
                profile = amp .* (par_pat(2)-1)/(pi*par_pat(1)^2) .* ...
                    (1+r.^2/par_pat(1)^2).^-par_pat(2) + ...
                    offset ; 
                amplitude = amp .* (par_pat(2)-1)/(pi*par_pat(1)^2) ;
            else
                profile = amp.*(1+r.^2/par_pat(1)^2).^-par_pat(2) + ...
                    offset ; 
                amplitude = amp ;
            end
        otherwise
            error([pattern_model, ...
                ': Unkown model for the pattern...']) ;
    end
    
    profile = OpDS*profile ;
end