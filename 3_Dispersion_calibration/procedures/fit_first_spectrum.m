%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to fit the first spectrum with the user supervision
% 
% Created: 03/07/2019 (mm/dd/yyyy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global list_amp_init
global list_y_int
global list_x_int

%% Initialization
% List of the positions
list_y_int = zeros(nb_spec_init, nb_lambda_cal) ;
list_x_int = zeros(nb_spec_init, nb_lambda_cal) ;

% List of the amplitudes
list_amp_init = val_max_wave*ones(nb_spec_init,nb_lambda_cal) ;

% Initialization of the pattern parameters
switch pattern_model.flag_profile
    case 'Gaussian'
        % One parameter: sigma
        par_pat = pattern_user_guess(1)*ones(1, nb_lambda_cal) ;
        if pattern_model.flag_norm
            list_amp_init = list_amp_init*2*pi*pattern_user_guess(1) ;
        end
    case 'Hexagon'
        % One parameter: sigma
        par_pat = pattern_user_guess(1)*ones(1, nb_lambda_cal) ;
        list_amp_init = list_amp_init*3^0.5*pattern_user_guess(1).^2 ;

    case 'Moffat'
        % Two parameters: alpha and beta
        par_pat = cat(3, ...
            pattern_user_guess(1)*ones(1, nb_lambda_cal), ...
            pattern_user_guess(2)*ones(1, nb_lambda_cal)) ; 
        if pattern_model.flag_norm
            list_amp_init = list_amp_init*pi* ...
                pattern_user_guess(1)^2/(pattern_user_guess(2)-1) ;
        end
end

% Scaling according to lambda if needed
if pattern_model.flag_scale
    for l = 1:nb_lambda_cal
        par_pat(1,l,1) = par_pat(1,l,1)*pattern_model.scale(l) ;
    end
end

%% Plot figures
fig_1 = figure(1) ;
subplot(1,2,1) ;
imshow(IFS_calib_wave, [val_min_wave, val_max_wave]) ;
plot_rectangle(fig_1, p_center, nb_center, 0, pix_IFS, ...
    [0.75, 0.75, 0.75], LineWidth) ;
subplot(1,2,2) ;
imshow(IFS_calib_spec, [val_min_spec, val_max_spec]) ;
plot_rectangle(fig_1, p_center, nb_center, 0, pix_IFS, ...
    [0.75, 0.75, 0.75], LineWidth) ;

fig_2 = figure(2) ;
subplot(1,2,1) ;
insert_calibration_pattern(fig_2, IFS_zoom_wave, pix_IFS, ...
    [], rad_ROI, pattern_model, par_pat, ...
    (list_amp_init-val_min_wave) / ...
    (val_max_wave-val_min_wave), 0, list_color) ;
title(['Select the center and a corner of a hexagon of the', ...
    ' first wavelength...']) ;
subplot(1,2,2) ;
insert_calibration_pattern(fig_2, IFS_zoom_spec, pix_IFS, ...
    [], rad_ROI, pattern_model, par_pat, ...
    (val_max_spec/val_max_wave*list_amp_init-val_min_spec) / ...
    (val_max_spec-val_min_spec), 0, list_color) ;
title(['Select the center and a corner of a hexagon of the', ...
    ' first wavelength...']) ;
disp(['Select the center and a corner of a hexagon of the', ...
    ' first wavelength... /!\ Insure that a complete hexagon', ...
    ' of each wavelength is in the field...']) ;

%% Select first 2 points
l = 1 ;
for hex = 1:2
    % User pointing
    [list_x_int(hex,l), list_y_int(hex,l)] = ...
        my_ginput(1, list_color(l,:)) ;
    list_y_int(hex,l) = list_y_int(hex,l)- floor(pix_zoom.nb_y/2+1) + ...
        p_center(1) ;
    list_x_int(hex,l) = list_x_int(hex,l)- floor(pix_zoom.nb_x/2+1) + ...
        p_center(2) ;

    % Loop on the calibrations
    for c = 1:nb_c
        % Refining position and amplitude
        for cc = 1:length(list_cal_pat_init)
            [list_y_int(hex,l), list_x_int(hex,l), ...
                list_amp_init(hex,l)] = ...
                find_pattern_pos_amp(IFS_calib_wave, ...
                [list_y_int(hex,l); list_x_int(hex,l); ...
                list_amp_init(hex,l)], ...
                rad_ROI, pix_IFS, option_opti_amp_pos, ...
                pattern_model, ...
                par_pat(1,l,:), 0, list_cal_pat_init{cc}, ...
                IFS_dark_wave, IFS_BP) ;
        end

        % Refining pattern parameter
        par_pat_update = get_pattern_param(par_pat, l) ;
        par_pat_update = find_pattern_par(IFS_calib_wave, ...
            par_pat_update, ...
            rad_ROI, pix_IFS, option_opti_par, pattern_model, ...
            list_y_int(hex,l), list_x_int(hex,l), ...
            list_amp_init(hex,l), 0, ...
            IFS_dark_wave, IFS_BP) ;
        par_pat = set_pattern_param(par_pat_update, par_pat, l) ;
    end

    % Display
    subplot(1,2,1) ;
    plot_rectangle(fig_2, ...
        [list_y_int(hex,l)-p_center(1), list_x_int(hex,l)-p_center(2)], ...
        [2*rad_ROI+1, 2*rad_ROI+1], 0, pix_zoom, list_color(l,:), ...
        LineWidth) ;

    subplot(1,2,2) ;
    plot_rectangle(fig_2, ...
        [list_y_int(hex,l)-p_center(1), list_x_int(hex,l)-p_center(2)], ...
        [2*rad_ROI+1, 2*rad_ROI+1], 0, pix_zoom, list_color(l,:), ...
        LineWidth) ;
end

%% Selecting the other wavelengths
hex = 1 ;
for l = 2:nb_lambda_cal
    % Display
    subplot(1,2,1) ;
    title(['Select the center the center of the hexagon', ...
        ' for the wavelength: ', num2str(l), '/', ...
        num2str(nb_lambda_cal)]) ;
    subplot(1,2,2) ;
    title(['Select the center the center of the hexagon', ...
        ' for the wavelength: ', num2str(l), '/', ...
        num2str(nb_lambda_cal)]) ;
    
    % User pointing
    [list_x_int(hex,l), list_y_int(hex,l)] = my_ginput(1, list_color(l,:)) ;
    list_y_int(hex,l) = list_y_int(hex,l) - floor(pix_zoom.nb_y/2+1) + ...
        p_center(1) ;
    list_x_int(hex,l) = list_x_int(hex,l) - floor(pix_zoom.nb_x/2+1) + ...
        p_center(2) ;

    % Loop on the calibrations
    for c = 1:nb_c
        % Refining position and amplitude
        for cc = 1:length(list_cal_pat_init)
            [list_y_int(hex,l), list_x_int(hex,l), ...
                list_amp_init(hex,l)] = ...
                find_pattern_pos_amp(IFS_calib_wave, ...
                [list_y_int(hex,l); list_x_int(hex,l); ...
                list_amp_init(hex,l)], ...
                rad_ROI, pix_IFS, option_opti_amp_pos, ...
                pattern_model, ...
                par_pat(1,l,:), 0, list_cal_pat_init{cc}, ...
                IFS_dark_wave, IFS_BP) ;
        end

        % Refining pattern parameter
        par_pat_update = get_pattern_param(par_pat, l) ;
        par_pat_update = find_pattern_par(IFS_calib_wave, ...
            par_pat_update, ...
            rad_ROI, pix_IFS, option_opti_par, pattern_model, ...
            list_y_int(hex,l), list_x_int(hex,l), ...
            list_amp_init(hex,l), 0, ...
            IFS_dark_wave, IFS_BP) ;
        par_pat = set_pattern_param(par_pat_update, par_pat, l) ;
    end

    % Display results
    for sub_plot_i = 1:2
        subplot(1,2,sub_plot_i) ;
        plot_rectangle(fig_2, ...
            [list_y_int(hex,l)-p_center(1), ...
            list_x_int(hex,l)-p_center(2)], ...
            [2*rad_ROI+1, 2*rad_ROI+1], 0, pix_zoom, list_color(l,:), ...
            LineWidth) ;
    end
    pause(0.1) ;
end

close(fig_1) ;
close(fig_2) ;

%% Estimation of the polynomial laws
% Estimation of the polynomial y-law
Pol = zeros(nb_lambda_cal, deg_pol_y+1) ;
for d = 0:deg_pol_y
    Pol(:,d+1) = list_lambda_cal.^d ;
end
coef_pol_y = (pinv(Pol)*list_y_int(1,:)')' ;


% Estimation of the polynomial x-law
Pol = zeros(nb_lambda_cal, deg_pol_x+1) ;
for d = 0:deg_pol_x
    Pol(:,d+1) = list_lambda_cal.^d ;
end
coef_pol_x = (pinv(Pol)*list_x_int(1,:)')' ;


% Estimation of the polynomial dif-law
Pol = zeros(nb_lambda_cal, deg_pol_dif+1) ;
for d = 0:deg_pol_dif
    Pol(:,d+1) = list_lambda_cal.^d ;
end
coef_pol_dif = (pinv(Pol)*par_pat')' ;


%% Extracting the region of interest
[list_i, list_j] = ...
    get_ROI(list_lambda, coef_pol_y', coef_pol_x', rad_ROI, pix_IFS) ;

[crop_y, crop_x] = get_extracted_pos(list_i, list_j, pix_IFS) ;

%% Forward models
% Selecting the wave calibration
Sel_wave = true(deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1 + ...
    nb_lambda + nb_lambda_cal, 1) ;
Sel_wave(deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1 + (1:nb_lambda)) = ...
    false ;
Sel_wave = LinOpSelector(Sel_wave) ;

% Selecting the spectrum calibration
Sel_spec = true(deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1 + ...
    nb_lambda + nb_lambda_cal, 1) ;
Sel_spec(deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1 + nb_lambda + ...
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
    'SPEC', 1, pattern_model)*Sel_spec ;

%% Building elementary cost
% Robust penalization options
option_robust_wave = [] ;
option_robust_wave.method = option_opti_spec_init.RL2_method ;
option_robust_wave.noise_model = option_opti_spec_init.noise_model ;
option_robust_wave.var_0 = option_opti_spec_init.var_0 ;
option_robust_wave.eta = option_opti_spec_init.eta ;
option_robust_wave.flag_s = IFS_W(list_i, list_j) .* ...
    option_opti_spec_init.flag_s ;
option_robust_spec = option_robust_wave ;

% Accounting for the dark current
option_robust_wave.var_0 = option_robust_wave.var_0 + ...
    option_robust_wave.eta*IFS_dark_wave(list_i, list_j) ;
option_robust_spec.var_0 = option_robust_spec.var_0 + ...
    option_robust_spec.eta*IFS_dark_spec(list_i, list_j) ;

% Building cost function
C_wave = CostRobustPenalization( ...
    Model_wave, IFS_calib_wave(list_i, list_j), ...
    option_robust_wave) ;
C_spec = CostRobustPenalization( ...
    Model_spec, IFS_calib_spec(list_i, list_j), ...
    option_robust_spec) ;
Sel_spec = false(deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1 + ...
    nb_lambda + nb_lambda_cal, 1) ;
Sel_spec(deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1 + (1:nb_lambda)) ...
    = true ;
Sel_spec = LinOpSelector(Sel_spec) ;
C_reg = CostL2([nb_lambda,1]) * LinOpGrad([nb_lambda,1], 1) * Sel_spec ;

% Normalization factors
norm_wave = 1/(length(list_i)*length(list_j)) ;
norm_spec = norm_wave ;
norm_lamb = 1/nb_lambda ;

%% Optimization
% Initial guess
par_spectrum = [ ...
    coef_pol_y'; ...
    coef_pol_x' ; ...
    coef_pol_dif'; ...
    2*val_max_spec*ones(nb_lambda,1); ...
    list_amp_init(1,:)'] ;

% Loop on the optimization
% Refining position and amplitude
for cc = 1:length(list_cal_spec_init)
    disp(['    Calibration of the first spectrum: ', num2str(cc), '/', ...
        num2str(length(list_cal_spec_init)), ...
        ' (', list_cal_spec_init{cc}, ')']) ;
    %% Initialization of the constraints
    cons_min = -Inf(size(par_spectrum)) ;
    cons_max = +Inf(size(par_spectrum)) ;
    
    % Constraints on the spectrum edge
    cons_max(deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1 + 1) = ...
        max_val_edge ;
    cons_max(deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1 + nb_lambda) = ...
        max_val_edge ;

    % Positivity constraints on the amplitude
    cons_min(deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1 + (1:nb_lambda)) ...
        = 0 ;
    cons_min(deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1 + nb_lambda + ...
        (1:nb_lambda_cal)) = 0 ;
    
    %% Building cost according the calibration type
    switch list_cal_spec_init{cc}
        case 'wave'
            % Constraints: fitting everything excepted the lamp spectrum
            cons_min(deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1 + ...
                (1:nb_lambda)) = par_spectrum(deg_pol_x+1 + ...
                deg_pol_y+1 + deg_pol_dif+1 + (1:nb_lambda)) ;
            cons_max(deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1 + ...
                (1:nb_lambda)) = par_spectrum(deg_pol_x+1 + ...
                deg_pol_y+1 + deg_pol_dif+1 + (1:nb_lambda)) ;
            
            % Global cost
            C = C_wave ;
            
        case 'spec'
            % Constraints: fitting only the lamp spectrum
            cons_min(1:(deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1)) = ...
                par_spectrum(1:(deg_pol_x+1 + deg_pol_y+1 + ...
                deg_pol_dif+1)) ;
            cons_min(deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1+ ...
                nb_lambda+(1:nb_lambda_cal)) = ...
                par_spectrum(deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1+ ...
                nb_lambda+(1:nb_lambda_cal)) ;
            cons_max(1:(deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1)) = ...
                par_spectrum(1:(deg_pol_x+1 + deg_pol_y+1 + ...
                deg_pol_dif+1)) ;
            cons_max(deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1+ ...
                nb_lambda+(1:nb_lambda_cal)) = ...
                par_spectrum(deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1+ ...
                nb_lambda+(1:nb_lambda_cal)) ;
            
            % Global cost
            C = norm_spec*C_spec ;
            if mu_L2_grad>0
                C = C + norm_lamb*mu_L2_grad / ...
                    (option_opti_spec_init.flag_s)^2*C_reg ;
            end
            
        case'both'
            % Global cost
            mu_wave_norm = mu_wave*(C_spec*par_spectrum) / ...
                (C_wave*par_spectrum) ;
            
            C = norm_wave*mu_wave_norm*C_wave + norm_spec*C_spec ;
            if mu_L2_grad>0
                C = C + norm_lamb*mu_L2_grad / ...
                    (option_opti_spec_init.flag_s)^2*C_reg ;
            end
    end
    
    %% Optimization
    par_spectrum_in = par_spectrum ;
    par_spectrum = run_Opti(C, cons_min, cons_max, ...
        option_opti_spec_init, par_spectrum) ;
    
%     %% Display results
%     dyn_res = 250 ;
%     % Wave
%     figure(42) ;
%     subplot(1,4,1) ;
%     imshow(IFS_calib_wave(list_i, list_j), [val_min_wave, val_max_wave]) ;
%     subplot(1,4,2) ;
%     imshow(Model_wave*par_spectrum, [val_min_wave, val_max_wave]) ;
%     subplot(1,4,3) ;
%     imshow(IFS_calib_wave(list_i, list_j)-Model_wave*par_spectrum, ...
%         [-dyn_res, dyn_res]) ;
%     subplot(1,4,4) ;
%     imshow(C_wave.computeW_(par_spectrum), [0, 1]) ;
%     
%     % Spec
%     figure(43) ;
%     subplot(1,4,1) ;
%     imshow(IFS_calib_spec(list_i, list_j), [val_min_spec, val_max_spec]) ;
%     subplot(1,4,2) ;
%     imshow(Model_spec*par_spectrum, [val_min_spec, val_max_spec]) ;
%     subplot(1,4,3) ;
%     imshow(IFS_calib_spec(list_i, list_j)-Model_spec*par_spectrum, ...
%         [-dyn_res, dyn_res]) ;
%     subplot(1,4,4) ;
%     imshow(C_spec.computeW_(par_spectrum), [0, 1]) ;
%     
%     % Lamp spectrum
%     figure(44) ;
%     hold on
%     plot(list_lambda, Sel_spec*par_spectrum) ;
%     axis([min(list_lambda), max(list_lambda), ...
%         0, max(Sel_spec*par_spectrum)]) ;
%     hold off
%     
%     disp([par_spectrum_in, par_spectrum]) ;
%     disp(par_spectrum_in-par_spectrum) ;
%     pause(0.1) ;
%     
% %     pause() ;
end

%% New parameters
coef_pol_y = par_spectrum(1:(deg_pol_y+1))' ;
coef_pol_x = par_spectrum(deg_pol_y+1 + (1:(deg_pol_x+1)))' ;
coef_pol_dif = par_spectrum(deg_pol_y+1 + deg_pol_x+1 + ...
    (1:(deg_pol_dif+1)))' ;
lamp_spec = par_spectrum(deg_pol_y+1 + deg_pol_x+1 + deg_pol_dif+1 + ...
    (1:nb_lambda)) ;
list_amp_init(1,:) = par_spectrum(deg_pol_y+1 + deg_pol_x+1 + ...
    deg_pol_dif+1 + nb_lambda + (1:nb_lambda_cal))' ;

mu_wave_norm = mu_wave*(C_spec*par_spectrum) / ...
    (C_wave*par_spectrum) ;

%% Display results
dyn_res = 250 ;
% Wave
fig_3 = figure(3) ;
subplot(2,8,1) ;
imshow(IFS_calib_wave(list_i, list_j), [val_min_wave, val_max_wave]) ;
subplot(2,8,2) ;
imshow(Model_wave*par_spectrum, [val_min_wave, val_max_wave]) ;
subplot(2,8,3) ;
imshow(IFS_calib_wave(list_i, list_j)-Model_wave*par_spectrum, ...
    [-dyn_res, dyn_res]) ;
subplot(2,8,4) ;
imshow(C_wave.computeW_(par_spectrum), [0, 1]) ;

% Spec
subplot(2,8,5) ;
imshow(IFS_calib_spec(list_i, list_j), [val_min_spec, val_max_spec]) ;
subplot(2,8,6) ;
imshow(Model_spec*par_spectrum, [val_min_spec, val_max_spec]) ;
subplot(2,8,7) ;
imshow(IFS_calib_spec(list_i, list_j)-Model_spec*par_spectrum, ...
    [-dyn_res, dyn_res]) ;
subplot(2,8,8) ;
imshow(C_spec.computeW_(par_spectrum), [0, 1]) ;

% Lamp spectrum
subplot(2,1,2) ;
plot(list_lambda, lamp_spec) ;
axis([min(list_lambda), max(list_lambda), 0, max(lamp_spec)]) ;
title('Is this guess correct? (y/n)') ;

str_user = input('Is this guess correct? (y/n)\n', 's') ;
y_user = strcmp(str_user, 'y') ;

close(fig_3) ;