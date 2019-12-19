%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to analyze the dark current with the shutter open
% 
% Created: 03/03/2019 (mm/dd/yyyy)
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Opening dark data
disp('Opening background darks...') ;
[nb_dark_BG, list_exp_time_dark_BG, list_avg_dark_BG, ...
    list_var_dark_BG] = ...
    open_list_files(data.path, data.dark_BG_name) ;
disp('Opening background darks: done!') ;


%% Fitting linear laws
[a_dark_BG, b_dark_BG] = ...
    fit_linear_law(list_exp_time_dark_BG, ...
    list_avg_dark_BG, 1./list_var_dark_BG, ...
    false) ;

[eta_dark_BG, sig_dark_BG] = ...
    fit_linear_law(list_exp_time_dark_BG, list_var_dark_BG, [], true) ;

%% Bad pixels identification
[BP_mask_dark_BG, BP_mask_dark_a_BG, BP_mask_dark_res_BG] = ...
    get_BP_mask(a_dark_BG, b_dark_BG, list_exp_time_dark_BG, ...
    list_avg_dark_BG, list_var_dark_BG, tol_MAD, size_med_filter) ;

%% Further cleaning
BP_mask_dark_BG(isnan(a_dark_BG)) = 0 ;
a_dark_BG(isnan(a_dark_BG)) = 1 ;
BP_mask_dark_BG(a_dark_BG<=0) = 0 ;
a_dark_BG(a_dark_BG<=0) = 0 ;

%% Saving
if save_lin_fits
    save_fits(a_dark_BG, 'a_dark_BG', [save_path, 'dark_BG/']) ;
    save_fits(b_dark_BG, 'b_dark_BG', [save_path, 'dark_BG/']) ;
    save_fits(eta_dark_BG, 'eta_dark_BG', [save_path, 'dark_BG/']) ;
    save_fits(sig_dark_BG, 'sig_dark_BG', [save_path, 'dark_BG/']) ;
    save_fits(uint8(BP_mask_dark_res_BG), 'BP_mask_dark_res_BG', ...
        [save_path, 'dark_BG/']) ;
    save_fits(uint8(BP_mask_dark_a_BG), 'BP_mask_dark_a_BG', ...
        [save_path, 'dark_BG/']) ;
end
save_fits(uint8(BP_mask_dark_BG), 'BP_mask_dark_BG', save_path) ;
% save_fits(uint8(BP_mask_dark_res_BG), 'BP_mask_dark_res_BG', save_path) ;
% save_fits(uint8(BP_mask_dark_a_BG), 'BP_mask_dark_a_BG', save_path) ;