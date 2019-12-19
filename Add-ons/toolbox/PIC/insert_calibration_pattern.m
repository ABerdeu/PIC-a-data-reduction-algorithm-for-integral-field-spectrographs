%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to insert the calibration pattern at the different wavelengths
% in the figure
% 
% Created: 09/03/2018 (mm/dd/yyyy)
% Modfied: 02/26/2019 (mm/dd/yyyy) Hexagon pattern
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #fig# figure in which the patterns must be inserted
%
% #pic# picture in which the patterns must be inserted
%
% #pix# pixel chracteristic
%   pix.nb_x -> number of the pixel on the x-axis
%   pix.nb_y -> number of the pixel on the y-axis
%   pix.dx -> size of the pixel on the x-axis
%   pix.dy -> size of the pixel on the y-axis
%
% #display_dyn# display dynamics (for imshow)
%
% #rad_cal# Number of pixel to plot the calibration patterns
%
% #pattern_model# flags to specify the pattern to fit
%   pattern_model.flag_norm -> Normalized?
%   pattern_model.flag_profile -> Gaussian / Moffat
%
% #par_pat# parameters of the pattern at the different wavelengths
%
% #amp_pat# amplitude of the pattern at the different wavelengths
%
% #off_pat# offset
%
% #flag_color# list of the colors of the different wavelengths
%
% #LineWidth# Line width
%
%%%%%%%%
% Ouput
%%%%%%%%
% #pic_out# picture with the inserted calibration patterns
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pic_out = insert_calibration_pattern(fig, pic, pix, dyn, ...
    rad_cal, pattern_model, par_pat, amp_pat, off_pat, flag_color, ...
    LineWidth)

    %% Options
    if nargin < 11 || isempty(LineWidth)
        LineWidth = 1 ;
    end

    %% Initialization
    nb_l = size(flag_color, 1) ;
    pic_out = pic ;
    nb_c = size(pic, 3) ;

    %% Indexes at which the calibration patterns are inserted
    ind_y_cal = 1+rad_cal+(-rad_cal:rad_cal) ;
    ind_x_cal = 1+rad_cal+(-rad_cal:rad_cal) ;
    
    %% Inserting patterns (with the parameters of the first pattern of the
    % list)
    for l = 1:nb_l
        par_pat_l = get_pattern_param(par_pat, l, 1) ;

        % Calibration pattern in the image corner
        cal_c = get_pattern_simulation( ...
                (-rad_cal:rad_cal)', (-rad_cal:rad_cal)', 0, 0, ...
                par_pat_l, ...
                median(amp_pat(:,l)), off_pat, pattern_model) ;
        for c = 1:nb_c
            pic_out(ind_y_cal+(l-1)*(2*rad_cal+1), ind_x_cal, c) = cal_c ;
        end
    end
    
    %% Display the picture with a rectange around the patterns
    imshow(pic_out, dyn) ;
    for l = 1:nb_l
        plot_rectangle(fig, ...
            [rad_cal+1+(l-1)*(2*rad_cal+1)-floor(pix.nb_y/2+1), ...
            rad_cal+1-floor(pix.nb_x/2+1)], ...
            [2*rad_cal+1, 2*rad_cal+1], 0, pix, flag_color(l,:), ...
            LineWidth) ;
    end
end