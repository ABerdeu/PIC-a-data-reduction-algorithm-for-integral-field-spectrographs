%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This macro is a procedure to automatically alternately refine the
% position and diffracion laws.
%
% Created: 03/15/2018 (mm/dd/yyyy)
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
% Modified: 10/14/2019 (mm/dd/yyyy) On-line deposit
% Author: Anthony Berdeu (Centre de Recherche en Astrophysique de Lyon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all ;
close all ;
clc ;

% Absolute path
abs_path = pwd ;
abs_path = [abs_path, '\'] ;

% Load functions and packages
restoredefaultpath

% Loading toolboxes
addpath(genpath('../Add-ons/')) ;

%% Parameters
path_cal = './results/2019_10_11_SAM/' ;

% Radius of the patterns in the projection matrix
rad_ROI_model = 5 ;

% Number of time that the refinement is performed
nb_cal = 10 ;


RL2_method = 'Cauchy' ;
noise_model = 'none' ;
flag_s = 350 ;

option_opti_dif.verbose = false ;
option_opti_dif.maxiter = 100 ;
option_opti_dif.method = 'fminsearch' ;

option_opti_lamp = option_opti_dif ;
option_opti_lamp.method = 'VMLMB' ;
option_opti_lamp.maxiter = 30 ;

option_opti_pos = option_opti_dif ;
option_opti_pos.method = 'fminsearch' ;
option_opti_pos.maxiter = 100 ;

mu_wave_pos = 50 ; % Weighting the wavelength positions compared to the
% spectral calibration in the cost function

nb_mu_wave_pos = 100 ; % Number of spectra on which the ratio between the
% costs of the wavelength and spectral calibration is estimated

% Format to save files
format_double = 'double' ;
format_uint8 = 'uint8' ;

%% Updating path
addpath(genpath(path_cal)) ;
if ~exist([path_cal, 'code/functions_refine/'], 'dir')
    save_path_aux = make_dir([path_cal, 'code/functions_refine/']) ;
end
copyfile('./refine_pos_dif_laws_parallel.m', [path_cal, 'code/']) ;
copyfile('./functions_refine', [path_cal, 'code/functions_refine/']) ;
addpath(genpath([path_cal, 'code/functions_refine/'])) ;

run parameters
save_path = [path_cal, 'refinement_pos_dif_par/'] ;
if ~exist(save_path, 'dir')
    save_path = make_dir(save_path) ;
end

% Auxiliary folders
save_path_coef = [save_path, 'coef/'] ;
if ~exist(save_path_coef, 'dir')
    save_path_coef = make_dir(save_path_coef) ;
end
save_path_wave = [save_path, 'wave/'] ;
if ~exist(save_path_wave, 'dir')
    save_path_wave = make_dir(save_path_wave) ;
end
save_path_spec = [save_path, 'spec/'] ;
if ~exist(save_path_spec, 'dir')
    save_path_spec = make_dir(save_path_spec) ;
end

%% Opening sensor calibration
% Mask on the bad pixels
IFS_BP = fitsread([data.sensor_calib, 'IFS_BP_mask.fits']) ;

% Sensor flat
IFS_sensor_flat = fitsread([data.sensor_calib, 'IFS_sensor_flat.fits']) ;
IFS_BP(isnan(IFS_sensor_flat)) = 0 ;
IFS_sensor_flat(isnan(IFS_sensor_flat)) = 1 ;
IFS_BP(IFS_sensor_flat<=0) = 0 ;
IFS_sensor_flat(IFS_sensor_flat<=0) = 1 ;


% Temporal dark current
IFS_dark_current = fitsread([data.sensor_calib, 'IFS_dark_BG.fits']) ;
IFS_BP(isnan(IFS_dark_current)) = 0 ;
IFS_dark_current(isnan(IFS_dark_current)) = 0 ;
IFS_BP(IFS_dark_current<=0) = 0 ;
IFS_dark_current(IFS_dark_current<=0) = 1 ;


% Mask on the scientific area
IFS_mask = fitsread([data.sensor_calib, 'IFS_mask.fits']) ;

% Variance law for the Poisson noise
if flag_MAD
    eta = load([data.sensor_calib, 'eta_MAD.txt'], '-ascii') ;
    var_0 = load([data.sensor_calib, 'sig_MAD.txt'], '-ascii') ;
else
    eta = load([data.sensor_calib, 'eta.txt'], '-ascii') ;
    var_0 = load([data.sensor_calib, 'sig.txt'], '-ascii') ;
end

option_opti_dif.var_0 = var_0 ;
option_opti_dif.eta = eta ;

% Opening data and correction by the sensor calibration
% Wave calibration file: correction by the dark current and the flat
IFS_calib_wave = fitsread(data.calib_wave) ;
IFS_calib_wave = mean(IFS_calib_wave,3) ; % Stacking the different frames
IFS_dark_wave = get_fits_exposure(data.calib_wave)*IFS_dark_current ./ ...
    IFS_sensor_flat ;
IFS_calib_wave = IFS_calib_wave./IFS_sensor_flat - IFS_dark_wave ;
IFS_BP(isnan(IFS_calib_wave)) = 0 ;
IFS_calib_wave(isnan(IFS_calib_wave)) = 0 ;
IFS_BP(IFS_calib_wave<=0) = 0 ;
IFS_calib_wave(IFS_calib_wave<=0) = 1 ;

% Wave calibration file: correction by the dark current and the flat
IFS_calib_spec = fitsread(data.calib_spec) ;
IFS_calib_spec = mean(IFS_calib_spec,3) ; % Stacking the different frames
IFS_dark_spec = get_fits_exposure(data.calib_spec)*IFS_dark_current ./ ...
    IFS_sensor_flat ;
IFS_calib_spec = IFS_calib_spec./IFS_sensor_flat - IFS_dark_spec ; 
IFS_BP(isnan(IFS_calib_spec)) = 0 ;
IFS_calib_spec(isnan(IFS_calib_spec)) = 0 ;
IFS_BP(IFS_calib_spec<=0) = 0 ;
IFS_calib_spec(IFS_calib_spec<=0) = 1 ;

% Pixel characteristic
[pix_IFS.nb_y, pix_IFS.nb_x, pix_IFS.nb_f] = size(IFS_calib_spec) ;
pix_IFS.dx = 1 ;
pix_IFS.dy = 1 ;

% Weigth on the defective pixels
IFS_W = IFS_BP ;
IFS_W(IFS_W==0) = Inf ;

% Display dynamics
val_min_wave = 0 ;
val_max_wave = IFS_calib_wave((IFS_BP>0)&(IFS_mask>0)) ;
val_max_wave = 10*median(val_max_wave(:)) ;
val_min_spec = 0 ;
val_max_spec = IFS_calib_spec((IFS_BP>0)&(IFS_mask>0)) ;
val_max_spec = 3*median(val_max_spec(:)) ;

%% Opening first calibration
coef_pol_y = load([path_cal, 'coef_pol_y.txt'], '-ascii') ;
coef_pol_x = load([path_cal, 'coef_pol_x.txt'], '-ascii') ;
coef_pol_dif = load([path_cal, 'coef_pol_dif.txt'], '-ascii') ;
trans = load([path_cal, 'trans.txt'], '-ascii') ;
list_amp = load([path_cal, 'list_amp.txt'], '-ascii') ;
list_lambda = load([path_cal, 'list_lambda.txt'], '-ascii') ;
list_lambda_cal = load([path_cal, 'list_lambda_cal.txt'], '-ascii') ;
lambda_0 = load([path_cal, 'lambda_0.txt'], '-ascii') ;
lamp_spec = load([path_cal, 'lamp_spec.txt'], '-ascii') ;

map_spec = fitsread([path_cal, 'map_spec.fits']) ;

% save_fits(IFS_calib_wave, 'IFS_calib_wave_corrected', save_path) ;
% save_fits(IFS_calib_spec, 'IFS_calib_spec_corrected', save_path) ;

save([save_path_coef, 'list_lambda.txt'], 'list_lambda', '-ascii', ...
    '-double') ;
save([save_path_coef, 'list_lambda_cal.txt'], ...
    'list_lambda_cal', '-ascii', '-double') ;
save([save_path_coef, 'lambda_0.txt'], 'lambda_0', '-ascii', '-double') ;

save([save_path_coef, 'coef_pol_y.txt'], 'coef_pol_y', '-ascii', ...
    '-double') ;
save([save_path_coef, 'coef_pol_x.txt'], 'coef_pol_x', '-ascii', ...
    '-double') ;
save([save_path_coef, 'coef_pol_dif.txt'], 'coef_pol_dif', '-ascii', ...
    '-double') ;
save([save_path_coef, 'trans.txt'], 'trans', '-ascii', '-double') ;
save([save_path_coef, 'list_amp.txt'], 'list_amp', '-ascii', '-double') ;
save([save_path_coef, 'lamp_spec.txt'], 'lamp_spec', '-ascii', '-double') ;
save_fits(map_spec, 'map_spec', save_path_coef) ;

% Relative list of lambda
list_lambda = list_lambda-lambda_0 ;
list_lambda_cal = list_lambda_cal-lambda_0 ;

nb_lambda = length(list_lambda) ;
nb_lambda_cal = length(list_lambda_cal) ;
list_color = jet(nb_lambda_cal) ;


%% Opening parallel jobs
% warning('off', 'MATLAB:imagesci:fitsinfo:fseekFailed') ;
% warning('off', 'MATLAB:mir_warning_maybe_uninitialized_temporary') ;

poolobj = gcp('nocreate');
if isempty(poolobj)
    try
        poolobj = parpool ;
        disp(['Connected to ', num2str(poolobj.NumWorkers), ...
        ' workers.']) ;
        nb_workers = poolobj.NumWorkers ;
    catch
        disp('Unable to start workers, license non available...') ;
        disp('No paralleling jobs...') ;
        nb_workers = 1 ;
    end
else
    disp(['Already connected to ', num2str(poolobj.NumWorkers), ...
        ' workers.']) ;
        nb_workers = poolobj.NumWorkers ;
end

% Path for the binary files
path_bin = [save_path, 'files_bin/'] ;
if ~exist(path_bin, 'dir')
    path_bin = make_dir(path_bin) ;
end

if nb_workers>1
    % Saving IFS_data and IFS_BP
    save_bin([path_bin, 'IFS_calib_wave.bin'], IFS_calib_wave, ...
        format_double) ;
    save_bin([path_bin, 'IFS_BP.bin'], IFS_BP, format_uint8) ;

    % Getting the pointers
    IFS_calib_wave_p = memmapfile([path_bin, 'IFS_calib_wave.bin'], ...
        'Format', format_double) ;
    IFS_BP_p = memmapfile([path_bin, 'IFS_BP.bin'], ...
        'Format', format_uint8) ;
else
    IFS_calib_wave_p = IFS_calib_wave ;
end


%% Hexagonal grid parameters
side_hex = ...
    (coef_pol_x(2:7,1) - coef_pol_x(1,1)) + ...
    (coef_pol_y(2:7,1) - coef_pol_y(1,1))*1i ;
side_hex = mean(side_hex,2) ;
side_hex = mean( abs(side_hex).* ...
    exp(1i*mod(angle(side_hex)*180/pi,60)*pi/180)) ;
theta_hex = mod(180/pi*angle(side_hex)+30,60) ;
pattern_model.theta = theta_hex ;
side_hex = abs(side_hex) ;
save([save_path_coef, '/theta_hex.txt'], 'theta_hex', '-ascii', ...
    '-double') ;
theta_hex = theta_hex-30 ;

%% Maximal number of hexagons and hexagon map
nb_rad = floor(side_hex/(2*rad_ROI)) ;
[pix_map.nb_y, pix_map.nb_x] = size(map_spec) ;
pix_map.dx = 1 ;
pix_map.dy = 1 ;
surf_hex = sum(IFS_mask(:)) ;
nb_spec = length(trans) ;

%% Display current results
% Estimation of the polynomial y-law
Pol_y_cal = zeros(nb_lambda_cal, deg_pol_y+1) ;
for d = 0:deg_pol_y
    Pol_y_cal(:,d+1) = list_lambda_cal.^d ;
end

% Estimation of the polynomial x-law
Pol_x_cal = zeros(nb_lambda_cal, deg_pol_x+1) ;
for d = 0:deg_pol_x
    Pol_x_cal(:,d+1) = list_lambda_cal.^d ;
end

% Estimation of the polynomial diffraction law
Pol_dif_cal = zeros(nb_lambda_cal, deg_pol_dif+1) ;
for d = 0:deg_pol_dif
    Pol_dif_cal(:,d+1) = list_lambda_cal.^d ;
end

% Loop on the spectra
list_y = zeros(nb_spec, nb_lambda_cal) ;
list_x = zeros(nb_spec, nb_lambda_cal) ;
for spec_aux = 1:nb_spec
    % Position of the patterns
    list_y(spec_aux,:) = Pol_y_cal*coef_pol_y(spec_aux,:)' ;
    list_x(spec_aux,:) = Pol_x_cal*coef_pol_x(spec_aux,:)' ;
end

fig_1 = figure(1) ;
subplot(1,2,1) ;
imshow(IFS_calib_wave, [val_min_wave, val_max_wave]) ;
for l = 1:nb_lambda_cal
    subplot(1,2,1) ;
    plot_dot(fig_1, ...
        [list_y(:,l), list_x(:,l)], ...
        pix_IFS, list_color(l,:), 10 * LineWidth) ;
end

subplot(1,2,2) ;
imshow(map_spec, []) ;
colormap(gca,jet) ;
title('Index map of the neighboring hexagons') ;

pause(0.1) ;

%% Loop on the calibrations
norm_spec = 1/(pix_IFS.nb_x*pix_IFS.nb_y) ;
norm_lamb = 1/nb_lambda ;
cons_min = zeros(nb_lambda, 1) ;
cons_max = Inf(nb_lambda, 1) ;

% Constraints on the spectrum edge
cons_max(1) = max_val_edge ;
cons_max(nb_lambda) = max_val_edge ;

option_robust_lamp = [] ;
option_robust_lamp.method = RL2_method ;
option_robust_lamp.noise_model = noise_model ;
option_robust_lamp.var_0 = var_0 ;
option_robust_lamp.eta = eta ;
option_robust_lamp.flag_s = IFS_W .* flag_s ;

% Initialization
coef_pol_y_ref = coef_pol_y ;
coef_pol_x_ref = coef_pol_x ;
coef_pol_dif_ref = coef_pol_dif ;
trans_ref = trans;
lamp_spec_ref = lamp_spec ;
list_amp_ref = list_amp ;

% Loop
for cc = 1:nb_cal
    disp(['Calibration: ', num2str(cc), '/', num2str(nb_cal)]) ;
    
    %% Building projection matrix
    spectra = [] ;
    spectra.list_lambda = list_lambda ;
    spectra.PSF_nb_x = 2*rad_ROI_model+1 ;
    spectra.PSF_nb_y = 2*rad_ROI_model+1 ;
    spectra.coef_pol_y = coef_pol_y_ref ;
    spectra.coef_pol_x = coef_pol_x_ref ;
    spectra.trans = trans_ref ;
    spectra.coef_pol_dif = coef_pol_dif_ref ;

    Spec_proj = LinOpSpecProjLocLaw(pix_IFS, spectra, pattern_model) ;
    
    %% Global residues
    on_lenslet = repmat(reshape(lamp_spec_ref, [1, nb_lambda]), ...
        [nb_spec, 1]) ;
    IFS_sim_spec = Spec_proj * on_lenslet ;
    save_fits(IFS_sim_spec, ['IFS_sim_spec_it_', num2str(cc-1)], ...
        save_path_spec) ;
    
    IFS_res_spec_p = IFS_calib_spec - IFS_sim_spec ;
    save_fits(IFS_res_spec_p, ...
        ['IFS_res_spec_it_', num2str(cc-1)], save_path_spec) ;
    
    %% Refining lamp_spec_res
    disp('    Refining lamp spectrum...') ;
    C_spec = CostRobustPenalization(Spec_proj * ...
        my_LinOpBroadcast([nb_spec, nb_lambda], 1), IFS_calib_spec, ...
        option_robust_lamp) ;

    % Regularization on the spectrum
    C_reg = CostL2([nb_lambda,1]) * LinOpGrad([nb_lambda,1], 1) ;

    % Global cost
    C = norm_spec*C_spec ;
    if mu_L2_grad>0
        C = C + norm_lamb*mu_L2_grad / ...
            (flag_s)^2*C_reg ;
    end

    % Optimization
    lamp_spec_ref_in = lamp_spec_ref ;
    lamp_spec_ref = run_Opti(C, cons_min, cons_max, ...
        option_opti_lamp, lamp_spec_ref_in) ;
    
    % Saving
    save([save_path_coef, 'lamp_spec_it_', num2str(cc), '.txt'], ...
        'lamp_spec_ref', '-ascii', '-double') ;
    
    on_lenslet = repmat(reshape(lamp_spec_ref, [1, nb_lambda]), ...
        [nb_spec, 1]) ;
    IFS_sim_spec = Spec_proj * on_lenslet ;
    IFS_res_spec_p = IFS_calib_spec - IFS_sim_spec ;
    if nb_workers>1
        % Saving IFS_res_spec_p
        save_bin([path_bin, 'IFS_res_spec.bin'], IFS_res_spec_p, ...
            format_double) ;

        % Getting the pointer
        IFS_res_spec_p = memmapfile([path_bin, 'IFS_res_spec.bin'], ...
            'Format', format_double) ;
    end
    
    %% Defining the ratio between the wavelength and the spectral
    % calibration
    if cc == 1
        disp(['    Initialization of the ratio between the ', ...
            'calibration files in the cost function...']) ;
        
        % Loop on a few spectrum to initialize the ratio
        list_cost = ones(nb_mu_wave_pos, 2) ;
        parfor spec = 1:nb_mu_wave_pos
            % Building the costs
            [C_wave, C_spec] = get_costs_wave_spec( ...
                coef_pol_y_ref(spec,:), coef_pol_x_ref(spec,:), ...
                coef_pol_dif_ref(spec,:), ...
                list_lambda_cal, list_lambda, trans_ref(spec), ...
                lamp_spec_ref, pattern_model, ...
                pix_IFS, rad_ROI, rad_ROI_model, ...
                IFS_calib_wave_p, IFS_res_spec_p, IFS_BP_p, ...
                RL2_method, noise_model, var_0, eta, flag_s) ;
            
            % Optimization
            par_in = [coef_pol_y_ref(spec,:)'; ...
                coef_pol_x_ref(spec,:)'; ...
                list_amp_ref(spec,:)'] ;

            % Update the variables
            list_cost(spec,:) = [C_wave*par_in, C_spec*par_in] ;
        end
        
        % Saving
        save([save_path_coef, 'list_cost_it_0.txt'], 'list_cost', ...
            '-ascii', '-double') ;
    end
    aux = list_cost(:,2)./list_cost(:,1) ;
    mu_wave_pos_norm = mu_wave_pos*median(aux(~isnan(aux))) ;
    save([save_path_coef, 'mu_wave_pos_norm_it_', num2str(cc), '.txt'], ...
        'mu_wave_pos_norm', '-ascii', '-double') ;
    
    %% Loop on the spectra
    [~, lastprint] = display_percentage('init', ...
        '    Fitting diffraction law') ;
    fprintf(repmat('\b', 1, lastprint));
    disp([': ', num2str(0, '%06.2f'), ' %']) ;
    lastprint = 11 ;
    save_fits(0, 'count', path_bin) ;
    delta_per = 0.1 ; % Percentage display
    delta_per = round(delta_per*nb_spec/100) ;
    list_cost = ones(nb_spec, 2) ;
    pctRunOnAll warning off
    parfor spec = 1:nb_spec
        %% Display the percentage
        display_percentage_par(path_bin, ...
            delta_per, nb_spec) ;
        
        %% Optimization of the position laws
        % Building the costs
        [C_wave, C_spec, Pol_y, Pol_x, Pol_dif, crop_y, crop_x, ...
            IFS_res_spec_ij, IFS_W_ij] = get_costs_wave_spec( ...
            coef_pol_y_ref(spec,:), coef_pol_x_ref(spec,:), ...
            coef_pol_dif_ref(spec,:), ...
            list_lambda_cal, list_lambda, trans_ref(spec), ...
            lamp_spec_ref, pattern_model, ...
            pix_IFS, rad_ROI, rad_ROI_model, ...
            IFS_calib_wave_p, IFS_res_spec_p, IFS_BP_p, ...
            RL2_method, noise_model, var_0, eta, flag_s) ;
        
        
        % Optimization
        par_in = [coef_pol_y_ref(spec,:)'; ...
            coef_pol_x_ref(spec,:)'; ...
            list_amp_ref(spec,:)'] ;
        C = mu_wave_pos_norm*C_wave + C_spec ;
        par_out = run_Opti(C, [], [], option_opti_pos, par_in) ;
        
        % Update the variables
        list_cost(spec,:) = [C_wave*par_out, C_spec*par_out] ;
        delta_y = ...
            par_out(1:(deg_pol_y+1))-coef_pol_y_ref(spec,:)' ;
        delta_y = delta_y(1) ;
        delta_x = ...
            par_out(deg_pol_y+1+(1:(deg_pol_x+1))) - ...
            coef_pol_x_ref(spec,:)' ;
        delta_x = delta_x(1) ;
        if (abs(delta_y) < rad_ROI) && (abs(delta_x) < rad_ROI)
            coef_pol_y_ref(spec,:) = par_out(1:(deg_pol_y+1)) ;
            coef_pol_x_ref(spec,:) = par_out(deg_pol_y+1 + ...
                (1:(deg_pol_x+1))) ;
            list_amp_ref(spec,:) = par_out(deg_pol_y+1+ deg_pol_x+1 + ...
                (1:nb_lambda_cal)) ;
        end
        
        
        %% Optimizing the diffraction law
        % Model of the lenslet
        % Local law of the lenslet
        p_y = Pol_y * coef_pol_y_ref(spec,:)' ;
        p_x = Pol_x * coef_pol_x_ref(spec,:)' ;

        % Elementary pattern
        Model_lambda = cell(nb_lambda, 1) ;
        for lambda = 1:nb_lambda
            Model_lambda{lambda} = get_pattern_operator( ...
                crop_y', crop_x', p_y(lambda), p_x(lambda), ...
                    [], [], 0, pattern_model) * LinOpMatrix( ...
                    [Pol_dif(lambda,:), 0 ;...
                    zeros(1, deg_pol_dif+1), lamp_spec_ref(lambda)]) ;
        end
        Model_lenslet = MapSummation(Model_lambda, ...
            ones(nb_lambda,1)) ;
        
        % Option of the robust penalization and cost
        option_robust = [] ;
        option_robust.method = RL2_method ;
        option_robust.noise_model = noise_model ;
        option_robust.var_0 = var_0 ;
        option_robust.eta = eta ;
        option_robust.flag_s = IFS_W_ij .* flag_s ;
        
        C = CostRobustPenalization( ...
            Model_lenslet, IFS_res_spec_ij, option_robust) ;
        
        % Optimization
        par_in = [coef_pol_dif_ref(spec,:)'; trans_ref(spec)] ;
        par_out = run_Opti(C, [], [], option_opti_dif, par_in) ;
        
        % Update the variables
        coef_pol_dif_ref(spec,:) = par_out(1:(deg_pol_dif+1)) ;
        trans_ref(spec) = par_out(deg_pol_dif+2) ;
        
    end
    pctRunOnAll warning on
    % Cleaning percentage
    display_percentage('exit', lastprint) ;
    
    %% Simulation of the wavelength calibration
    % Loop on the spectra
    list_y = zeros(nb_spec, nb_lambda_cal) ;
    list_x = zeros(nb_spec, nb_lambda_cal) ;
    list_dif = zeros(nb_spec, nb_lambda_cal) ;
    for spec_aux = 1:nb_spec
        % Position of the patterns
        list_y(spec_aux,:) = Pol_y_cal*coef_pol_y_ref(spec_aux,:)' ;
        list_x(spec_aux,:) = Pol_x_cal*coef_pol_x_ref(spec_aux,:)' ;
        
        % Diffraction law of the pattern
        list_dif(spec_aux,:) = Pol_dif_cal*coef_pol_dif_ref(spec_aux,:)' ;
    end
    
    IFS_sim_wave = simulate_IFS(pix_IFS, pattern_model, list_y, ...
        list_x, list_dif, list_amp_ref, 0, 2*rad_ROI, true) ;
    
    %% Saving
    save([save_path_coef, 'list_amp_it_', num2str(cc), '.txt'], ...
        'list_amp_ref', '-ascii', '-double') ;
    save([save_path_coef, 'coef_pol_y_it_', num2str(cc), '.txt'], ...
        'coef_pol_y_ref', '-ascii', '-double') ;
    save([save_path_coef, 'coef_pol_x_it_', num2str(cc), '.txt'], ...
        'coef_pol_x_ref', '-ascii', '-double') ;
    save([save_path_coef, 'coef_pol_dif_it_', num2str(cc), '.txt'], ...
        'coef_pol_dif_ref', '-ascii', '-double') ;
    save([save_path_coef, 'trans_it_', num2str(cc), '.txt'], ...
        'trans_ref', '-ascii', '-double') ;
    save([save_path_coef, 'list_cost_it_', num2str(cc), '.txt'], ...
        'list_cost', '-ascii', '-double') ;
    save_fits(IFS_sim_wave, ['IFS_sim_wave_it_', num2str(cc)], ...
        save_path_wave) ;
    save_fits(IFS_calib_wave - IFS_sim_wave, ...
        ['IFS_res_wave_it_', num2str(cc)], save_path_wave) ;
end


%% Last global simulation
% Building projection matrix
spectra = [] ;
spectra.list_lambda = list_lambda ;
spectra.PSF_nb_x = 2*rad_ROI_model+1 ;
spectra.PSF_nb_y = 2*rad_ROI_model+1 ;
spectra.coef_pol_y = coef_pol_y_ref ;
spectra.coef_pol_x = coef_pol_x_ref ;
spectra.trans = trans_ref ;
spectra.coef_pol_dif = coef_pol_dif_ref ;

Spec_proj = LinOpSpecProjLocLaw(pix_IFS, spectra, pattern_model) ;


% Global residues
on_lenslet = repmat(reshape(lamp_spec_ref, [1, nb_lambda]), ...
    [nb_spec, 1]) ;
IFS_sim_spec = Spec_proj * on_lenslet ;
save_fits(IFS_sim_spec, ['IFS_sim_spec_it_', num2str(nb_cal)], ...
    save_path_spec) ;
save_fits(IFS_calib_spec - IFS_sim_spec, ...
    ['IFS_res_spec_it_', num2str(nb_cal)], save_path_spec) ;

%% Removing auxiliary parallel files
remove_dir(path_bin) ;