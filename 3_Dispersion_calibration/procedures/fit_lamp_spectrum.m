%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to fit the lamp spectrum using neighboring spectra
% 
% Created: 03/08/2019 (mm/dd/yyyy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global list_y_int
global list_x_int
global list_amp_init

%% Initialization
% Coefficients of the law
coef_pol_y_init = repmat(coef_pol_y, [nb_spec_init, 1]) ;
coef_pol_x_init = repmat(coef_pol_x, [nb_spec_init, 1]) ;
coef_pol_dif_init = repmat(coef_pol_dif, [nb_spec_init, 1]) ;
list_amp_init = repmat(list_amp_init(1,:), [nb_spec_init, 1]) ;
list_cost = ones(nb_spec_init, 2) ;
list_cost(1,2) = mu_wave_norm ;

% List of the transmission
trans_init = ones(nb_spec_init, 1) ;

%% Hexagonal grid characteristic
side_hex = (list_x_int(2,1)-list_x_int(1,1)) + ...
    (list_y_int(2,1)-list_y_int(1,1))*1i ;
theta_hex = mod(180/pi*angle(side_hex)+30,60)-30 ;
side_hex = abs(side_hex) ;

% Orientation of the lenslet -> +30° compared to the hexagonal grid
pattern_model.theta = theta_hex+30 ;

%% Defining the hexagone map
% Radius of a hexagon number in the hexagon map
nb_rad = floor(side_hex/(2*rad_ROI)) ;
map_hex = zeros(pix_IFS.nb_y+2*round(side_hex), ...
    pix_IFS.nb_x+2*round(side_hex), ...
    'uint16') ;
[pix_map.nb_y, pix_map.nb_x] = size(map_hex) ;
pix_map.dx = 1 ;
pix_map.dy = 1 ;
ind_hex = pos2ind([coef_pol_y_init(1,1), coef_pol_x_init(1,1)], pix_map) ;
map_hex(ind_hex(1,1)+(-nb_rad*rad_ROI:nb_rad*rad_ROI), ...
    ind_hex(1,2)+(-nb_rad*rad_ROI:nb_rad*rad_ROI)) = 1 ;

%% Loop to find the neighboring spectra
spec_current = 1 ;
spec = 2 ;
while spec<nb_spec_init+1
    % Automatic rough pointing
    % List of the local hexagon corners
    list_neigh = get_neigh_hex(map_hex, ...
        [coef_pol_y_init(spec_current,1), ...
        coef_pol_x_init(spec_current,1)], side_hex, ...
        theta_hex) ;

    % Neighboring corners
    list_c = (1:6)' ;
    list_c = list_c(list_neigh>0) ;
    list_c_new = setdiff((1:6)', list_c) ;

    if ~isempty(list_c_new)
        disp(['    Fitting spectrum: ', num2str(spec), '/', ...
            num2str(nb_spec_init)]) ;
        %% Creating the local spectrum
        c_new = list_c_new(1) ;
        [x_c, y_c] = rot_2D(theta_hex+(c_new-1)*60, ...
            side_hex, 0) ;
        
        % Updating position
        coef_pol_x_init(spec,1) = x_c + coef_pol_x_init(spec_current,1) ;
        coef_pol_y_init(spec,1) = y_c + coef_pol_y_init(spec_current,1) ;
        
        %% Fitting spectrum
        [coef_pol_y_init(spec,:), coef_pol_x_init(spec,:), ...
            coef_pol_dif_init(spec,:), trans_init(spec), ...
            list_amp_init(spec,:), ~, ~, ~, ~,...
            list_cost(spec,1), list_cost(spec,2)...
            ] = fit_calib_spectrum( ...
            coef_pol_y_init(spec,:), coef_pol_x_init(spec,:), ...
            coef_pol_dif_init(spec,:), trans_init(spec), ...
            list_amp_init(spec,:), list_lambda, list_lambda_cal, ...
            rad_ROI, pix_IFS, lamp_spec, ...
            pattern_model, IFS_calib_wave, IFS_calib_spec, ...
            IFS_dark_wave, IFS_dark_spec, ...
            IFS_W, mu_wave_norm, option_opti_spec_init) ;
        
        %% Update the index map
        ind_hex = pos2ind( ...
            [coef_pol_y_init(spec,1), coef_pol_x_init(spec,1)], ...
            pix_map) ;
        map_hex(...
            ind_hex(1)+(-nb_rad*rad_ROI:nb_rad*rad_ROI), ...
            ind_hex(2)+(-nb_rad*rad_ROI:nb_rad*rad_ROI)) = ...
            spec ;
        spec = spec+1 ;
    else
        spec_current = spec_current+1 ;
    end
end

%% New normalization for the cost
mu_wave_norm = mu_wave*median(list_cost(:,2)./list_cost(:,1)) ;

%% Extracting position of the calibration patterns
min_i = [] ;
max_i = [] ;
min_j = [] ;
max_j = [] ;

% Estimation of the polynomial y-law
Pol_y = zeros(nb_lambda_cal, deg_pol_y+1) ;
for d = 0:deg_pol_y
    Pol_y(:,d+1) = list_lambda_cal.^d ;
end

% Estimation of the polynomial x-law
Pol_x = zeros(nb_lambda_cal, deg_pol_x+1) ;
for d = 0:deg_pol_x
    Pol_x(:,d+1) = list_lambda_cal.^d ;
end

% Loop on the spectra
for spec = 1:nb_spec_init
    % Position of the patterns
    list_y_int(spec,:) = Pol_y*coef_pol_y_init(spec,:)' ;
    list_x_int(spec,:) = Pol_x*coef_pol_x_init(spec,:)' ;
    
    % Extracting the region of interest
    [list_i, list_j] = ...
        get_ROI(list_lambda, coef_pol_y_init(spec,:)', ...
        coef_pol_x_init(spec,:)', rad_ROI, pix_IFS) ;
    
    % Update of the extraction parameters
    min_i = min([min_i, list_i]) ;
    max_i = max([max_i, list_i]) ;
    min_j = min([min_j, list_j]) ;
    max_j = max([max_j, list_j]) ;
end

%% Display fitted positions
fig_1 = figure(1) ;
subplot(1,2,1) ;
imshow(IFS_calib_wave, [val_min_wave, val_max_wave]) ;
plot_rectangle(fig_1, p_center, nb_center, 0, pix_IFS, ...
    [0.75, 0.75, 0.75]) ;
subplot(1,2,2) ;
imshow(IFS_calib_spec, [val_min_spec, val_max_spec]) ;
plot_rectangle(fig_1, p_center, nb_center, 0, pix_IFS, ...
    [0.75, 0.75, 0.75]) ;

fig_2 = figure(2) ;
subplot(1,2,1) ;
insert_calibration_pattern(fig_2, IFS_zoom_wave, pix_IFS, ...
    [], rad_ROI, pattern_model, par_pat, ...
    (list_amp_init(1,:)-val_min_wave) / ...
    (val_max_wave-val_min_wave), 0, list_color) ;
subplot(1,2,2) ;
insert_calibration_pattern(fig_2, IFS_zoom_spec, pix_IFS, ...
    [], rad_ROI, pattern_model, par_pat, ...
    (val_max_spec/val_max_wave*list_amp_init(1,:)-val_min_spec) / ...
    (val_max_spec-val_min_spec), 0, list_color) ;


subplot(1,2,2) ;
plot_circle(fig_2, ...
    [coef_pol_y_init(:,1)-p_center(1), ...
    coef_pol_x_init(:,1)-p_center(2)], ...
    rad_ROI, pix_zoom, [0.5, 1, 0.5], LineWidth) ;

for l = 1:nb_lambda_cal
    subplot(1,2,1) ;
    plot_circle(fig_2, ...
        [list_y_int(:,l)-p_center(1), list_x_int(:,l)-p_center(2)], ...
        rad_ROI, pix_zoom, list_color(l,:), LineWidth) ;
end

pause(0.1) ;

%% Extracting the region of interest
list_i = min_i:max_i ;
list_j = min_j:max_j ;
[crop_y, crop_x] = get_extracted_pos(list_i, list_j, pix_IFS) ;

%% Initialization
% Robust penalization options
option_robust_wave = [] ;
option_robust_wave.method = option_opti_spec.RL2_method ;
option_robust_wave.noise_model = option_opti_spec.noise_model ;
option_robust_wave.var_0 = option_opti_spec.var_0 ;
option_robust_wave.eta = option_opti_spec.eta ;
option_robust_wave.flag_s = IFS_W(list_i, list_j) .* ...
    option_opti_spec.flag_s ;
option_robust_spec = option_robust_wave ;

% Accounting for the dark current
option_robust_wave.var_0 = option_robust_wave.var_0 + ...
    option_robust_wave.eta*IFS_dark_wave(list_i, list_j) ;
option_robust_spec.var_0 = option_robust_spec.var_0 + ...
    option_robust_spec.eta*IFS_dark_spec(list_i, list_j) ;

%% Loop on the optimization
% Refining position and amplitude
for cc = 1:length(list_cal_spec_esti)
    disp(['    Joint calibration of the first spectra: ', ...
        num2str(cc), '/', num2str(length(list_cal_spec_esti)), ...
        ' (', list_cal_spec_esti{cc}, ')']) ;
    
    %% Building cost according the calibration type
    switch list_cal_spec_esti{cc}
        case 'wave'
            %% Wave
            coef_pol_y_init_aux = coef_pol_y_init' ;
            coef_pol_x_init_aux = coef_pol_x_init' ;
            coef_pol_dif_init_aux = coef_pol_dif_init' ;
            list_amp_init_aux = list_amp_init' ;

            % Initial guess
            par_spectrum = [ ...
                coef_pol_y_init_aux(:); ...
                coef_pol_x_init_aux(:) ; ...
                coef_pol_dif_init_aux(:); ...
                list_amp_init_aux(:)] ;
            
            % Wave calibration
            [Model_wave] = get_model_spec(crop_y', crop_x', ...
                deg_pol_y, deg_pol_x, deg_pol_dif, list_lambda_cal, ...
                nb_spec_init, 'SPEC', trans_init, pattern_model) ;
            
            %% Optimization
            C_wave = CostRobustPenalization( ...
                Model_wave, IFS_calib_wave(list_i, list_j), ...
                option_robust_wave) ;
            
            % Optimization
            par_spectrum_in = par_spectrum ;
            par_spectrum = run_Opti(C_wave, [], [], ...
                option_opti_spec, ...
                par_spectrum_in) ;
            
            %% Update the variables
            % coef_pol_y
            coef_pol_y_init = par_spectrum( ...
                1:nb_spec_init*(deg_pol_y+1)) ;
            coef_pol_y_init = reshape(coef_pol_y_init, ...
                [deg_pol_y+1, nb_spec_init])' ;
            
            % coef_pol_x
            coef_pol_x_init = par_spectrum( ...
                nb_spec_init*(deg_pol_x+1) + ...
                (1:nb_spec_init*(deg_pol_x+1))) ;
            coef_pol_x_init = reshape(coef_pol_x_init, ...
                [deg_pol_x+1, nb_spec_init])' ;
            
            % coef_pol_dif
            coef_pol_dif_init = par_spectrum( ...
                nb_spec_init*(deg_pol_y+1 + deg_pol_x+1) + ...
                (1:nb_spec_init*(deg_pol_dif+1))) ;
            coef_pol_dif_init = reshape(coef_pol_dif_init, ...
                [deg_pol_dif+1, nb_spec_init])' ;
            
            % list_amp_init
            list_amp_init = par_spectrum( ...
                nb_spec_init * ...
                (deg_pol_y+1 + deg_pol_x+1 + deg_pol_dif+1) + ...
                (1:nb_spec_init*nb_lambda_cal)) ;
            list_amp_init = reshape(list_amp_init, ...
                [nb_lambda_cal, nb_spec_init])' ;
            
%             %% Display results
%             dyn_res = 700 ;
%             % Wave
%             figure(42) ;
%             subplot(1,4,1) ;
%             imshow(IFS_calib_wave(list_i, list_j), ...
%                 [val_min_wave, val_max_wave]) ;
%             subplot(1,4,2) ;
%             imshow(Model_wave*par_spectrum, [val_min_wave, val_max_wave]) ;
%             subplot(1,4,3) ;
%             imshow(IFS_calib_wave(list_i, list_j)-...
%                 Model_wave*par_spectrum, ...
%                 [-dyn_res, dyn_res]) ;
%             subplot(1,4,4) ;
%             imshow(C_wave.computeW_(par_spectrum), [0, 1]) ;
%             pause() ;

            
        case 'spec'
            %% Spec
            % Building projection matrix
            sensor = [] ;
            sensor.dx = 1 ;
            sensor.dy = 1 ;
            sensor.nb_y = length(list_i) ;
            sensor.nb_x = length(list_j) ;
            
            spectra = [] ;
            spectra.trans = trans_init ;
            spectra.list_lambda = list_lambda ;
            spectra.PSF_nb_x = 5*rad_ROI ;
            spectra.PSF_nb_y = 5*rad_ROI ;
            spectra.coef_pol_y = coef_pol_y_init ;
            spectra.coef_pol_y(:,1) = spectra.coef_pol_y(:,1) - ...
                crop_y(floor(sensor.nb_y/2+1)) ;
            spectra.coef_pol_x = coef_pol_x_init ;
            spectra.coef_pol_x(:,1) = spectra.coef_pol_x(:,1) - ...
                crop_x(floor(sensor.nb_x/2+1)) ;
            spectra.coef_pol_dif = coef_pol_dif_init ;
            
            Model_spec = LinOpSpecProjLocLaw(sensor, ...
                spectra, pattern_model) * ...
                my_LinOpBroadcast([nb_spec_init, nb_lambda], 1) ;

            % Initial guess
            par_spectrum = lamp_spec ;
            
            % Initialization of the constraints
            cons_min = zeros(size(lamp_spec)) ;
            cons_max = Inf(size(lamp_spec)) ;
            
            % Constraints on the spectrum edge
            cons_max(1) = max_val_edge ;
            cons_max(nb_lambda) = max_val_edge ;
            
            %% Optimization
            C_spec = CostRobustPenalization( ...
                Model_spec, IFS_calib_spec(list_i, list_j), ...
                option_robust_spec) ;
            
            % Regularization on the spectrum
            C_reg = CostL2([nb_lambda,1]) * ...
                LinOpGrad([nb_lambda,1], 1) ;
            
            % Normalization factors
            norm_spec = 1/(length(list_i)*length(list_j)) ;
            norm_lamb = 1/nb_lambda ;
            
            % Global cost
            C = norm_spec*C_spec ;
            if mu_L2_grad>0
                C = C + norm_lamb*mu_L2_grad / ...
                    (option_opti_spec.flag_s)^2*C_reg ;
            end
            
            % Optimization
            par_spectrum_in = par_spectrum ;
            par_spectrum = run_Opti(C, cons_min, cons_max, ...
                option_opti_spec, par_spectrum_in) ;
            
            %% Update the variables
            % coef_pol_dif
            lamp_spec = par_spectrum ;
            
%             %% Display results
%             par_plot = par_spectrum ;
%             
%             dyn_res = 700 ;
%             % Wave
%             figure(43) ;
%             subplot(1,4,1) ;
%             imshow(IFS_calib_spec(list_i, list_j), ...
%                 [val_min_spec, val_max_spec]) ;
%             subplot(1,4,2) ;
%             imshow(Model_spec*par_plot, ...
%                 [val_min_spec, val_max_spec]) ;
%             subplot(1,4,3) ;
%             imshow(IFS_calib_spec(list_i, list_j) - ...
%                 Model_spec*par_plot, ...
%                 [-dyn_res, dyn_res]) ;
%             subplot(1,4,4) ;
%             imshow(C_spec.computeW_(par_plot), [0, 1]) ;
% 
%             
%             % Lamp spectrum
%             figure(44) ;
%             hold on
%             plot(list_lambda, par_plot) ;
%             axis([min(list_lambda), max(list_lambda), ...
%                 0, max(par_plot)]) ;
%             hold off
% 
%             disp([par_spectrum_in, par_spectrum]) ;
%             disp(par_spectrum_in-par_spectrum) ;
%             pause(0.1) ;
%             pause() ;
            
        case'both'
            %% both
            coef_pol_y_init_aux = coef_pol_y_init' ;
            coef_pol_x_init_aux = coef_pol_x_init' ;
            coef_pol_dif_init_aux = coef_pol_dif_init' ;
            list_amp_init_aux = list_amp_init' ;

            % Initial guess
            par_spectrum = [ ...
                coef_pol_y_init_aux(:); ...
                coef_pol_x_init_aux(:) ; ...
                coef_pol_dif_init_aux(:); ...
                list_amp_init_aux(:); ...
                trans_init] ;            
            
            % Wave calibration
            Sel_wave = true(nb_spec_init * ...
                (deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1 + ...
                nb_lambda_cal + 1), 1) ;
            Sel_wave(nb_spec_init * ...
                (deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1 + ...
                nb_lambda_cal)+(1:nb_spec_init)) = ...
                false ;
            Sel_wave = LinOpSelector(Sel_wave) ;
            [Model_wave] = get_model_spec(crop_y', crop_x', ...
                deg_pol_y, deg_pol_x, deg_pol_dif, list_lambda_cal, ...
                nb_spec_init, 'SPEC', trans_init, pattern_model)*Sel_wave ;
            
            % Spectral calibration
            Sel_spec = true(nb_spec_init * ...
                (deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1 + ...
                nb_lambda_cal + 1), 1) ;
            Sel_spec(nb_spec_init * ...
                (deg_pol_x+1 + deg_pol_y+1 + deg_pol_dif+1) + ...
                (1:nb_spec_init*nb_lambda_cal)) = ...
                false ;
            Sel_spec = LinOpSelector(Sel_spec) ;
            [Model_spec] = get_model_spec(crop_y', crop_x', ...
                deg_pol_y, deg_pol_x, deg_pol_dif, list_lambda, ...
                nb_spec_init, 'TRANS', lamp_spec, pattern_model)*Sel_spec ;
            
            %% Optimization
            C_wave = CostRobustPenalization( ...
                Model_wave, IFS_calib_wave(list_i, list_j), ...
                option_robust_wave) ;
            
            C_spec = CostRobustPenalization( ...
                Model_spec, IFS_calib_spec(list_i, list_j), ...
                option_robust_spec) ;
            
            % Normalization factors
            norm_wave = 1/(length(list_i)*length(list_j)) ;
            norm_spec = norm_wave ;
            
            % Global cost
            C = norm_wave*mu_wave_norm*C_wave + norm_spec*C_spec ;
            
            % Optimization
            par_spectrum_in = par_spectrum ;
            par_spectrum = run_Opti(C, [], [], ...
                option_opti_spec, par_spectrum_in) ;
            
            %% Update the variables
            % coef_pol_y
            coef_pol_y_init = par_spectrum( ...
                1:nb_spec_init*(deg_pol_y+1)) ;
            coef_pol_y_init = reshape(coef_pol_y_init, ...
                [deg_pol_y+1, nb_spec_init])' ;
            
            % coef_pol_x
            coef_pol_x_init = par_spectrum( ...
                nb_spec_init*(deg_pol_x+1) + ...
                (1:nb_spec_init*(deg_pol_x+1))) ;
            coef_pol_x_init = reshape(coef_pol_x_init, ...
                [deg_pol_x+1, nb_spec_init])' ;
            
            % coef_pol_dif
            coef_pol_dif_init = par_spectrum( ...
                nb_spec_init*(deg_pol_y+1 + deg_pol_x+1) + ...
                (1:nb_spec_init*(deg_pol_dif+1))) ;
            coef_pol_dif_init = reshape(coef_pol_dif_init, ...
                [deg_pol_dif+1, nb_spec_init])' ;
            
            % list_amp_init
            list_amp_init = par_spectrum( ...
                nb_spec_init * ...
                (deg_pol_y+1 + deg_pol_x+1 + deg_pol_dif+1) + ...
                (1:nb_spec_init*nb_lambda_cal)) ;
            list_amp_init = reshape(list_amp_init, ...
                [nb_lambda_cal, nb_spec_init])' ;
            
            % trans_init
            trans_init = par_spectrum( ...
                nb_spec_init * ...
                (deg_pol_y+1 + deg_pol_x+1 + deg_pol_dif+1 + ...
                nb_lambda_cal) + (1:nb_spec_init)) ;
            
%             %% Display results
%             par_plot = par_spectrum ;
%             dyn_res = 700 ;
%             % Wave
%             figure(42) ;
%             subplot(1,4,1) ;
%             imshow(IFS_calib_wave(list_i, list_j), ...
%                 [val_min_wave, val_max_wave]) ;
%             subplot(1,4,2) ;
%             imshow(Model_wave*par_plot, [val_min_wave, val_max_wave]) ;
%             subplot(1,4,3) ;
%             imshow(IFS_calib_wave(list_i, list_j)-Model_wave*par_plot, ...
%                 [-dyn_res, dyn_res]) ;
%             subplot(1,4,4) ;
%             imshow(C_wave.computeW_(par_plot), [0, 1]) ;
%             
%             % Spec
%             figure(43) ;
%             subplot(1,4,1) ;
%             imshow(IFS_calib_spec(list_i, list_j), ...
%                 [val_min_spec, val_max_spec]) ;
%             subplot(1,4,2) ;
%             imshow(Model_spec*par_plot, [val_min_spec, val_max_spec]) ;
%             subplot(1,4,3) ;
%             imshow(IFS_calib_spec(list_i, list_j)-Model_spec*par_plot, ...
%                 [-dyn_res, dyn_res]) ;
%             subplot(1,4,4) ;
%             imshow(C_spec.computeW_(par_plot), [0, 1]) ;
%             
%             pause() ;
    end
end

save([save_path, 'lamp_spec.txt'], 'lamp_spec', '-ascii', '-double') ;

%% Display results
dyn_res = 750 ;

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

% save_fits(IFS_calib_spec(list_i, list_j), 'CALIB', save_path) ;
% save_fits(Model_spec*par_spectrum, 'SIMU', save_path) ;
% save_fits(IFS_calib_spec(list_i, list_j) - Model_spec*par_spectrum, ...
%     'RES', save_path) ;

% Lamp spectrum
subplot(2,1,2) ;
plot(list_lambda, lamp_spec) ;
axis([min(list_lambda), max(list_lambda), 0, max(lamp_spec)]) ;

pause(0.1) ;
