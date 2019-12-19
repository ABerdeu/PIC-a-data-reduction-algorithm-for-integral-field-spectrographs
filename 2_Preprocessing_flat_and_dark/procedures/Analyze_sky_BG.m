%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to analyze the sky background
% 
% Created: 03/04/2019 (mm/dd/yyyy)
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Opening dark data
disp('Opening sky backgrounds...') ;
[nb_sky_BG, list_exp_time_sky_BG, list_avg_sky_BG, ...
    list_var_sky_BG] = ...
    open_list_files(data.path, data.sky_BG_name) ;
disp('Opening sky backgrounds: done!') ;


%% Fitting linear laws
[a_sky_BG, b_sky_BG] = ...
    fit_linear_law(list_exp_time_sky_BG, ...
    list_avg_sky_BG, [], ...
    false) ;


%% Further cleaning
BP_mask(isnan(a_sky_BG)) = 0 ;
a_sky_BG(isnan(a_sky_BG)) = 1 ;
BP_mask(a_sky_BG<=0) = 0 ;
a_sky_BG(a_sky_BG<=0) = 0 ;


%% Saving
if save_lin_fits
    save_fits(a_sky_BG, 'a_sky_BG', [save_path, 'sky_BG/']) ;
    save_fits(b_sky_BG, 'b_sky_BG', [save_path, 'sky_BG/']) ;
end