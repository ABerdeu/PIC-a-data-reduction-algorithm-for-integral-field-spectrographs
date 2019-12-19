%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to analyze the sensor flat
% 
% Created: 03/03/2019 (mm/dd/yyyy)
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Opening dark data
disp('Opening flats...') ;
[nb_flat, list_exp_time_flat, list_avg_flat, list_var_flat] = ...
    open_list_files(data.path, data.flat_name) ;
disp('Opening flats: done!') ;

%% Correction of the dark current and indexing with the median exposure
list_med_flat = zeros(nb_flat,1) ;
for flat = 1:nb_flat
    list_avg_flat(:,:,flat) = list_avg_flat(:,:,flat) - ...
        (list_exp_time_flat(flat).*a_dark+b_dark) ;
    med_flat_aux = list_avg_flat(:,:,flat) ;
    list_med_flat(flat) = median(med_flat_aux(:)) ;
end
clear('med_flat_aux') ;

%% Fitting linear laws
[a_flat, b_flat] = ...
    fit_linear_law(list_med_flat, list_avg_flat, 1./list_var_flat, ...
    false) ;

[eta_flat, sig_flat] = ...
    fit_linear_law(list_med_flat, list_var_flat, [], false) ;

% Removing NaN
a_flat(isnan(a_flat)) = 0 ;
b_flat(isnan(b_flat)) = 0 ;
eta_flat(isnan(eta_flat)) = 0 ;
sig_flat(isnan(sig_flat)) = 0 ;


%% Bad pixels identification
[BP_mask_flat, BP_mask_flat_a, BP_mask_flat_res] = ...
    get_BP_mask(a_flat, b_flat, list_med_flat, ...
    list_avg_flat, list_var_flat, tol_MAD, size_med_filter) ;


%% Further cleaning
BP_mask_flat(isnan(a_flat)) = 0 ;
a_flat(isnan(a_flat)) = 1 ;
BP_mask_flat(a_flat<=0) = 0 ;
a_flat(a_flat<=0) = 1 ;


%% Saving
if save_lin_fits
    save_fits(a_flat, 'a_flat', [save_path, 'flat/']) ;
    save_fits(b_flat, 'b_flat', [save_path, 'flat/']) ;
    save_fits(eta_flat, 'eta_flat', [save_path, 'flat/']) ;
    save_fits(sig_flat, 'sig_flat', [save_path, 'flat/']) ;
    save_fits(uint8(BP_mask_dark_res_BG), 'BP_mask_dark_res_BG', ...
        [save_path, 'flat/']) ;
    save_fits(uint8(BP_mask_flat_a), 'BP_mask_flat_a', ...
        [save_path, 'flat/']) ;
end
save_fits(uint8(BP_mask_flat), 'BP_mask_flat', save_path) ;
% save_fits(uint8(BP_mask_flat_res), 'BP_mask_flat_res', save_path) ;
% save_fits(uint8(BP_mask_flat_a), 'BP_mask_flat_a', save_path) ;