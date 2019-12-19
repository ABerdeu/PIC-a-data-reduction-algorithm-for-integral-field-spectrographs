%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to fit a spectrum for the wavelength and spectral calibration
% files
% 
% Created: 03/11/2019 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #[coef_pol_y_in, coef_pol_x_in, coef_pol_dif_in]# initial value for the
% polynomial laws
% 
% #trans_in# initial value for the transmission
%
% #list_amp_in# initial fitted amplitude for the wavelength calibration
%
% #list_lambda# list of the simulation wavelength
%
% #list_lambda_cal# list of the calibration wavelength
%
% #rad_ROI# radius on which the patterns are fitted
%
% #pix# pixel caracteristics
%
% #lamp_spec# spectrum of the lamp
%
% #pattern_model# individual pattern model
%
% #W# Weight matrix to exclude defective pixels
%
% #[Calib_wave, Calib_spec]# Calibration files for the wavelength and the
% lamp spectrum
%
% #[Dark_wave, Dark_spec]# Dark currents for the files
%
% #mu_wave# Hyper-parameters to weight the influence of the wavelength
% calibration compared to the spectral calibration.
%   -> mu_wave = 0, only the diffraction law is fitted
%
% #option_opti_m# option for the opimization
%   -> RL2_method
%   -> noise_model
%   -> (var_0)
%   -> (eta)
%   -> flag_s
%   -> structure containing the options for the optimizer
%
%%%%%%%%
% Ouput
%%%%%%%%
% #[coef_pol_y_out, coef_pol_x_out, coef_pol_dif_out]# final value for the
% polynomial laws
%
% #trans_out# final value for the transmission
%
% #list_amp_out# final fitted amplitude for the wavelength calibration
%
% #[list_i, list_j]# positions in the grid on which the spectrum is fitted
%
% #[Sim_wave, Sim_spec]# Simulated fit on the region of interest
%
% #[c_wave, c_spec]# The final value of the costs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [coef_pol_y_out, coef_pol_x_out, coef_pol_dif_out, ...
    trans_out, list_amp_out, list_i, list_j, Sim_wave, Sim_spec, ...
    c_wave, c_spec] = ...
    fit_calib_spectrum( ...
    coef_pol_y_in, coef_pol_x_in, coef_pol_dif_in, trans_in, ...
    list_amp_in, list_lambda, list_lambda_cal, rad_ROI, pix, lamp_spec, ...
    pattern_model, Calib_wave, Calib_spec, Dark_wave, Dark_spec, ...
    W, mu_wave, option_opti_m)

    %% Initialization
    deg_pol_y = length(coef_pol_y_in) - 1 ;
    deg_pol_x = length(coef_pol_x_in) - 1 ;
    deg_pol_dif = length(coef_pol_dif_in) - 1 ;
    nb_lambda_cal = length(list_lambda_cal) ;
    
    %% Extracting the region of interest
    [list_i, list_j] = ...
        get_ROI(list_lambda, coef_pol_y_in', ...
        coef_pol_x_in', rad_ROI, pix) ;

    [crop_y, crop_x] = get_extracted_pos(list_i, list_j, pix) ;

    %% Forward models
    % Selecting the wave calibration
    Sel_wave = true(deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1 + ...
        1 + nb_lambda_cal, 1) ;
    Sel_wave(deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1 + 1) = ...
        false ;
    Sel_wave = LinOpSelector(Sel_wave) ;

    % Selecting the spectrum calibration
    Sel_spec = true(deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1 + ...
        1 + nb_lambda_cal, 1) ;
    Sel_spec(deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1 + 1 + ...
        (1:nb_lambda_cal)) = ...
        false ;
    Sel_spec = LinOpSelector(Sel_spec) ;

    % Wave calibration
    [Model_wave] = get_model_spec(crop_y', crop_x', ...
        deg_pol_y, deg_pol_x, deg_pol_dif, list_lambda_cal, 1, ...
        'SPEC', 1, pattern_model)*Sel_wave ;

    % Spectrum calibration
    [Model_spec] = get_model_spec(crop_y', crop_x', ...
        deg_pol_y, deg_pol_x, deg_pol_dif, list_lambda, 1, ...
        'TRANS', lamp_spec, pattern_model)*Sel_spec ;


    %% Building elementary cost
    % Robust penalization options
    option_robust_wave = [] ;
    option_robust_wave.method = option_opti_m.RL2_method ;
    option_robust_wave.noise_model = ...
        option_opti_m.noise_model ;
    option_robust_wave.var_0 = option_opti_m.var_0 ;
    option_robust_wave.eta = option_opti_m.eta ;
    option_robust_wave.flag_s = W(list_i, list_j) .* ...
        option_opti_m.flag_s ;
    option_robust_spec = option_robust_wave ;

    % Accounting for the dark current
    option_robust_wave.var_0 = option_robust_wave.var_0 + ...
        option_robust_wave.eta*Dark_wave(list_i, list_j) ;
    option_robust_spec.var_0 = option_robust_spec.var_0 + ...
        option_robust_spec.eta*Dark_spec(list_i, list_j) ;

    % Building cost function
    C_wave = CostRobustPenalization( ...
        Model_wave, Calib_wave(list_i, list_j), ...
        option_robust_wave) ;
    C_spec = CostRobustPenalization( ...
        Model_spec, Calib_spec(list_i, list_j), ...
        option_robust_spec) ;

    %% Optimization
    % Initial guess
    par_spectrum = [ ...
        coef_pol_y_in'; ...
        coef_pol_x_in' ; ...
        coef_pol_dif_in'; ...
        trans_in; ...
        list_amp_in'] ;

    %% Initialization of the constraints and global cost
    cons_min = -Inf(size(par_spectrum)) ;
    cons_max = +Inf(size(par_spectrum)) ;

    % Limitation on the divergence of the positions
    cons_min(1) = coef_pol_y_in(1)-rad_ROI ;
    cons_min(deg_pol_y+1 + 1) = coef_pol_x_in(1)-rad_ROI ;
    cons_max(1) = coef_pol_y_in(1)+rad_ROI ;
    cons_max(deg_pol_y+1 + 1) = coef_pol_x_in(1)+rad_ROI ;

    % Cost
    if mu_wave==0
        % Fitting only the diffraction law
        
        % Fixing the position law
        cons_min(1:(deg_pol_y+1 + deg_pol_x+1)) = ...
            par_spectrum(1:(deg_pol_y+1 + deg_pol_x+1)) ;
        cons_max(1:(deg_pol_y+1 + deg_pol_x+1)) = ...
            par_spectrum(1:(deg_pol_y+1 + deg_pol_x+1)) ;
        
        % Fixing the amplitudes
        cons_min( deg_pol_y+1 + deg_pol_x+1 + deg_pol_dif+1 + ...
            1 + (1:nb_lambda_cal)) = ...
            par_spectrum(deg_pol_y+1 + deg_pol_x+1 + deg_pol_dif+1 + ...
            1 + (1:nb_lambda_cal)) ;
        cons_max( deg_pol_y+1 + deg_pol_x+1 + deg_pol_dif+1 + ...
            1 + (1:nb_lambda_cal)) = ...
            par_spectrum(deg_pol_y+1 + deg_pol_x+1 + deg_pol_dif+1 + ...
            1 + (1:nb_lambda_cal)) ;
        
        % Cost on the spectral calibration
        C = C_spec ;
    else
%%%
%   -> mu_wave > 0, the diffraction law is fixed
%         % Fixing the diffraction law
%         cons_min( deg_pol_y+1 + deg_pol_x+1 + ...
%             (1:(deg_pol_dif+1))) = ...
%             par_spectrum( deg_pol_y+1 + deg_pol_x+1 + ...
%             (1:(deg_pol_dif+1))) ;
%         cons_max( deg_pol_y+1 + deg_pol_x+1 + ...
%             (1:(deg_pol_dif+1))) = ...
%             par_spectrum( deg_pol_y+1 + deg_pol_x+1 + ...
%             (1:(deg_pol_dif+1))) ;
%%%
        C = mu_wave*C_wave + C_spec ;
    end
    

% %     %% Display initialization
% %     dyn_res = 750 ;
% %     % Wave
% %     figure(40) ;
% %     subplot(1,4,1) ;
% %     imshow(Calib_wave(list_i, list_j), [0, 2000]) ;
% %     subplot(1,4,2) ;
% %     imshow(Model_wave*par_spectrum, [0, 2000]) ;
% %     subplot(1,4,3) ;
% %     imshow(Calib_wave(list_i, list_j)-Model_wave*par_spectrum, ...
% %         [-dyn_res, dyn_res]) ;
% %     subplot(1,4,4) ;
% %     imshow(C_wave.computeW_(par_spectrum), [0, 1]) ;
% %     
% %     % Spec
% %     figure(41) ;
% %     subplot(1,4,1) ;
% %     imshow(Calib_spec(list_i, list_j), [0, 9000]) ;
% %     subplot(1,4,2) ;
% %     imshow(Model_spec*par_spectrum, [0, 9000]) ;
% %     subplot(1,4,3) ;
% %     imshow(Calib_spec(list_i, list_j)-Model_spec*par_spectrum, ...
% %         [-dyn_res, dyn_res]) ;
% %     subplot(1,4,4) ;
% %     imshow(C_spec.computeW_(par_spectrum), [0, 1]) ;
% %     pause(0.1) ;
% %     

    %% Optimization
    par_spectrum = run_Opti(C, cons_min, cons_max, ...
        option_opti_m, par_spectrum) ;

    %% New parameters
    coef_pol_y_out = par_spectrum(1:(deg_pol_y+1))' ;
    coef_pol_x_out = par_spectrum(deg_pol_y+1 + ...
        (1:(deg_pol_x+1)))' ;
    coef_pol_dif_out = par_spectrum(deg_pol_y+1 + ...
        deg_pol_x+1 + (1:(deg_pol_dif+1)))' ;
    trans_out = par_spectrum(deg_pol_y+1 + deg_pol_x+1 + ...
        deg_pol_dif+1 + 1) ;
    list_amp_out = par_spectrum(deg_pol_y+1 + ...
        deg_pol_x+1 + deg_pol_dif+1 + 1 + (1:nb_lambda_cal))' ;
    
    % Simulation
    Sim_wave = Model_wave*par_spectrum ;
    Sim_spec = Model_spec*par_spectrum ;
    
    c_wave = C_wave*par_spectrum ;
    c_spec = C_spec*par_spectrum ;
    
% %     %% Display results
% %     dyn_res = 750 ;
% %     % Wave
% %     figure(42) ;
% %     subplot(1,4,1) ;
% %     imshow(Calib_wave(list_i, list_j), [0, 2000]) ;
% %     subplot(1,4,2) ;
% %     imshow(Model_wave*par_spectrum, [0, 2000]) ;
% %     subplot(1,4,3) ;
% %     imshow(Calib_wave(list_i, list_j)-Model_wave*par_spectrum, ...
% %         [-dyn_res, dyn_res]) ;
% %     subplot(1,4,4) ;
% %     imshow(C_wave.computeW_(par_spectrum), [0, 1]) ;
% %     
% %     % Spec
% %     figure(43) ;
% %     subplot(1,4,1) ;
% %     imshow(Calib_spec(list_i, list_j), [0, 9000]) ;
% %     subplot(1,4,2) ;
% %     imshow(Model_spec*par_spectrum, [0, 9000]) ;
% %     subplot(1,4,3) ;
% %     imshow(Calib_spec(list_i, list_j)-Model_spec*par_spectrum, ...
% %         [-dyn_res, dyn_res]) ;
% %     subplot(1,4,4) ;
% %     imshow(C_spec.computeW_(par_spectrum), [0, 1]) ;
% %     
% %     pause(0.1) ;
end

