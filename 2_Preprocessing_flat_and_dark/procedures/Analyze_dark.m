%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to analyze the dark current
% 
% Created: 02/25/2019 (mm/dd/yyyy)
% Modified: 03/03/2019 (mm/dd/yyyy)
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Opening dark data
disp('Opening darks...') ;
[nb_dark, list_exp_time_dark, list_avg_dark, list_var_dark] = ...
    open_list_files(data.path, data.dark_name) ;
disp('Opening darks: done!') ;


%% Fitting linear laws
[a_dark, b_dark] = ...
    fit_linear_law(list_exp_time_dark, list_avg_dark, 1./list_var_dark, ...
    false) ;

[eta_dark, sig_dark] = ...
    fit_linear_law(list_exp_time_dark, list_var_dark, [], true) ;


%% Bad pixels identification
[BP_mask_dark, BP_mask_dark_a, BP_mask_dark_res] = ...
    get_BP_mask(a_dark, b_dark, list_exp_time_dark, list_avg_dark, ...
    list_var_dark, tol_MAD, size_med_filter) ;

%% Further cleaning
BP_mask_dark(isnan(a_dark)) = 0 ;
a_dark(isnan(a_dark)) = 1 ;
BP_mask_dark(a_dark<=0) = 0 ;
a_dark(a_dark<=0) = 0 ;


%% Saving
if save_lin_fits
    save_fits(a_dark, 'a_dark', [save_path, 'dark/']) ;
    save_fits(b_dark, 'b_dark', [save_path, 'dark/']) ;
    save_fits(eta_dark, 'eta_dark', [save_path, 'dark/']) ;
    save_fits(sig_dark, 'sig_dark', [save_path, 'dark/']) ;
    save_fits(uint8(BP_mask_dark_res), 'BP_mask_dark_res', ...
        [save_path, 'dark/']) ;
    save_fits(uint8(BP_mask_dark_a), 'BP_mask_dark_a', ...
        [save_path, 'dark/']) ;
end
save_fits(uint8(BP_mask_dark), 'BP_mask_dark', save_path) ;
% save_fits(uint8(BP_mask_dark_res), 'BP_mask_dark_res', save_path) ;
% save_fits(uint8(BP_mask_dark_a), 'BP_mask_dark_a', save_path) ;