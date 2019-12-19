%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to estimate the position shift between the two files
% 
% Created: 11/20/2018 (mm/dd/yyyy)
% Created: 03/13/2019 (mm/dd/yyyy) Adaptation to local polynomial laws and
% autocalibration on the data
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Estimation of the position shift...') ;

save_path_shift = make_dir([save_path, 'estimate_shift/']) ;
save_path_shift_it = make_dir([save_path_shift, 'iterations/']) ;

%% Loading variables
% Mask on the bad pixels
IFS_BP = fitsread([path_calib_sensor, 'IFS_BP_mask.fits']) ;

% Sensor flat
IFS_sensor_flat = fitsread([path_calib_sensor, 'IFS_sensor_flat.fits']) ;

% Temporal dark current
IFS_dark_current = fitsread([path_calib_sensor, 'IFS_dark_BG.fits']) ;

% Mask on the scientific area
IFS_mask = fitsread([path_calib_sensor, 'IFS_mask.fits']) ;

% Variance law for the Poisson noise
if flag_MAD
    eta = load([path_calib_sensor, 'eta_MAD.txt'], '-ascii') ;
    var_0 = load([path_calib_sensor, 'sig_MAD.txt'], '-ascii') ;
else
    eta = load([path_calib_sensor, 'eta.txt'], '-ascii') ;
    var_0 = load([path_calib_sensor, 'sig.txt'], '-ascii') ;
end

% Data
IFS_data = fitsread(path_data) ;
IFS_data = mean(IFS_data,3) ; % Stacking the different frames
% % save_fits(IFS_data, 'IFS_data_not_corrected', save_path_shift) ;
IFS_dark_wave = exp_time*IFS_dark_current ./ ...
    IFS_sensor_flat ;
IFS_data = IFS_data./IFS_sensor_flat - IFS_dark_wave ;
IFS_BP(isnan(IFS_data)) = 0 ;
IFS_data(isnan(IFS_data)) = 0 ;
IFS_BP(IFS_data<=0) = 0 ;
IFS_data(IFS_data<=0) = 1 ;

save_fits(IFS_data, 'IFS_data_corrected', save_path_shift) ;

% Spectra map
spec_map = fitsread([path_calib_spec, 'map_spec.fits']) ;

% Weigth on the defective pixels
IFS_W = IFS_BP ;
IFS_W(IFS_W==0) = Inf ;

% Pixel characteristic
[pix_IFS.nb_y, pix_IFS.nb_x, pix_IFS.nb_f] = size(IFS_data) ;
pix_IFS.dx = 1 ;
pix_IFS.dy = 1 ;

% Pixel characteristic
[pix_map.nb_y, pix_map.nb_x, pix_map.nb_f] = size(spec_map) ;
pix_map.dx = 1 ;
pix_map.dy = 1 ;

% Plotting values
val_min = 0 ;
val_max = IFS_data((IFS_BP>0)&(IFS_mask>0)) ;
val_max = 10*median(val_max(:)) ;

%% Selecting regions of interest
% Position of the regions of interest
list_x_ROI = ((1:nb_ROI)-(nb_ROI+1)/2)/(nb_ROI+1) * ...
    min(pix_IFS.nb_y, pix_IFS.nb_x) ;
[list_x_ROI, list_y_ROI] = meshgrid(list_x_ROI, list_x_ROI) ;
list_x_ROI = list_x_ROI(:) ;
list_y_ROI = list_y_ROI(:) ;
[list_x_ROI, list_y_ROI] = rot_2D(theta_square, list_x_ROI, list_y_ROI) ;
sel_ROI = list_x_ROI>-0.95*pix_IFS.nb_x/2 & ...
    list_x_ROI<0.95*pix_IFS.nb_x/2 & ...
    list_y_ROI>-0.95*pix_IFS.nb_y/2 & ...
    list_y_ROI<0.95*pix_IFS.nb_y/2 ;
list_x_ROI = list_x_ROI(sel_ROI) ;
list_y_ROI = list_y_ROI(sel_ROI) ;
nb_tot_ROI = length(list_x_ROI) ;

% Display
fig_1 = figure(1) ;
set(fig_1,'rend','painters','pos', [fig_pos, fig_size]) ;
imshow(IFS_data, [val_min, val_max]) ;
plot_rectangle(fig_1, [list_y_ROI, list_x_ROI], ...
    [2*rad_ROI_shift+1, 2*rad_ROI_shift+1], 0, pix_IFS, ...
    [1, 0.5, 0.5]) ;
title('Selected regions of interest') ;
set(gca,'FontSize',FontSize_axis);
title('Selected regions of interest', ...
    'FontSize', FontSize) ;
set(findall(gca,'-property','Interpreter'),'Interpreter','Latex') ;
set(findall(gcf,'-property','Interpreter'),'Interpreter','Latex') ;
pause(1) ;
saveas(fig_1, [save_path_shift, 'Selected regions of interest'], 'png') ;


%% Extracting regions of interest
% Finding lenslet of each region of interest
list_ind = pos2ind([list_y_ROI, list_x_ROI], pix_map) ;
list_lenslet = cell(nb_tot_ROI, 1) ;
list_nb_lenslet = zeros(nb_tot_ROI, 1) ;
for roi = 1:nb_tot_ROI
    list_ind_roi = spec_map( ...
        list_ind(roi,1)+(-rad_ROI_shift:rad_ROI_shift), ...
        list_ind(roi,2)+(-rad_ROI_shift:rad_ROI_shift)) ;
    list_lenslet{roi} = setdiff(list_ind_roi(:), 0) ;
    list_nb_lenslet(roi) = length(list_lenslet{roi}) ;
end

% Extracting regions of interest
list_ind = pos2ind([list_y_ROI, list_x_ROI], pix_IFS) ;
pic_ROI = zeros(2*rad_ROI_shift+1, 2*rad_ROI_shift+1, nb_tot_ROI) ;
IFS_W_ROI = pic_ROI ;
coef_pol_y_ROI = cell(nb_tot_ROI, 1) ;
coef_pol_x_ROI = cell(nb_tot_ROI, 1) ;
coef_pol_dif_ROI = cell(nb_tot_ROI, 1) ;
trans_ROI = cell(nb_tot_ROI, 1) ;
for roi = 1:nb_tot_ROI
    % Extracting ROI
    pic_ROI(:,:,roi) = IFS_data( ...
        list_ind(roi,1)+(-rad_ROI_shift:rad_ROI_shift), ...
        list_ind(roi,2)+(-rad_ROI_shift:rad_ROI_shift)) ;
    IFS_W_ROI(:,:,roi) = IFS_W( ...
        list_ind(roi,1)+(-rad_ROI_shift:rad_ROI_shift), ...
        list_ind(roi,2)+(-rad_ROI_shift:rad_ROI_shift)) ;

    % Extraction coefficients
    coef_pol_y_ROI{roi} = coef_pol_y(list_lenslet{roi},:) ;
    coef_pol_x_ROI{roi} = coef_pol_x(list_lenslet{roi},:) ;
    coef_pol_dif_ROI{roi} = coef_pol_dif(list_lenslet{roi},:) ;
    trans_ROI{roi} = trans(list_lenslet{roi},:) ;

    % Centering the coefficients
    coef_pol_y_ROI{roi}(:,1) = coef_pol_y_ROI{roi}(:,1) + ...
        floor(pix_IFS.nb_y/2)+1 - list_ind(roi,1) ;
    coef_pol_x_ROI{roi}(:,1) = coef_pol_x_ROI{roi}(:,1)  + ...
        floor(pix_IFS.nb_x/2)+1 - list_ind(roi,2) ;
end

%% Building the equivalent data
pic_data_ROI = zeros(nb_tot_ROI^0.5*(2*rad_ROI_shift+1) + ...
    nb_tot_ROI^0.5 +1 ) ;
for roi = 1:nb_tot_ROI
    j_roi = 1+mod(roi-1, nb_tot_ROI^0.5) ;
    i_roi = 1+floor(mod((roi-1)/nb_tot_ROI^0.5, nb_tot_ROI^0.5)) ;
    
    j_roi = j_roi + (j_roi-1)*(2*rad_ROI_shift+1) + ...
        (1:(2*rad_ROI_shift+1)) ;
    i_roi = i_roi + (i_roi-1)*(2*rad_ROI_shift+1) + ...
        (1:(2*rad_ROI_shift+1)) ;
    
    pic_data_ROI(i_roi,j_roi) = pic_ROI(:,:,roi) ;
end
save_fits(pic_data_ROI, 'pic_data_ROI', save_path_shift) ;

%% Loop on the calibrations
% Initialization
delta_pos = [0;0] ;
spec_roi = val_max*ones(nb_lambda, nb_tot_ROI) ;
sensor_ROI = [] ;
sensor_ROI.dx = 1 ;
sensor_ROI.dy = 1 ;
sensor_ROI.nb_y = 2*rad_ROI_shift+1 ;
sensor_ROI.nb_x = 2*rad_ROI_shift+1 ;
pic_ROI_spec = pic_ROI ;
pic_W_coef_spec = pic_ROI ;

% Normalization factors
norm_spec = 1/(2*rad_ROI_shift+1)^2 ;
norm_lamb = 1/nb_lambda ;
cons_min = zeros(nb_lambda, 1) ;
cons_max = Inf(nb_lambda, 1) ;
% Constraints on the spectrum edge
cons_max(1) = max_val_edge ;
cons_max(nb_lambda) = max_val_edge ;

% Positions of the pixels on the sensor
sensor_x = (-rad_ROI_shift:rad_ROI_shift)' ;
sensor_y = sensor_x ;

% Polynomial laws
% y-law polynomial
deg_pol_y = size(coef_pol_y, 2)-1 ;
Pol_y = zeros(nb_lambda, deg_pol_y+1) ;
for d = 0:deg_pol_y
    Pol_y(:,d+1) = (list_lambda-lambda_0).^d ;
end

% x-law polynomial
deg_pol_x = size(coef_pol_x, 2)-1 ;
Pol_x = zeros(nb_lambda, deg_pol_x+1) ;
for d = 0:deg_pol_x
    Pol_x(:,d+1) = (list_lambda-lambda_0).^d ;
end

% Diffraction polynomial
deg_pol_dif = size(coef_pol_dif, 2)-1 ;
Pol_dif = zeros(nb_lambda, deg_pol_dif+1) ;
for d = 0:deg_pol_dif
    Pol_dif(:,d+1) = (list_lambda-lambda_0).^d ;
end

% Loop on the calibrations
for cc = 1:nb_cal
    disp(['Calibration: ', num2str(cc), '/', num2str(nb_cal)]) ;
    %% Spectra fitting
    % Loop on the regions of interest
    pic_sim_ROI = zeros(nb_tot_ROI^0.5*(2*rad_ROI_shift+1) + ...
        nb_tot_ROI^0.5 +1 ) ;
    pic_W_coef_ROI = zeros(nb_tot_ROI^0.5*(2*rad_ROI_shift+1) + ...
        nb_tot_ROI^0.5 +1 ) ;
    for roi = 1:nb_tot_ROI
        disp(['   Spectra fitting ', num2str(roi), '/', ...
            num2str(nb_tot_ROI)]) ;
        % Option of the robust penalization
        option_robust = [] ;
        option_robust.method = option_opti_spec.RL2_method ;
        option_robust.noise_model = option_opti_spec.noise_model ;
        option_robust.var_0 = var_0 ;
        option_robust.eta = eta ;
        option_robust.flag_s = IFS_W_ROI(:,:,roi) .* ...
            option_opti_spec.flag_s ;

        % Building projection matrix
        spectra_ROI = spectra ;
        spectra_ROI.trans = trans_ROI{roi} ;
        spectra_ROI.coef_pol_y = coef_pol_y_ROI{roi} ;
        spectra_ROI.coef_pol_y(:,1) = spectra_ROI.coef_pol_y(:,1) + ...
            +delta_pos(1) ;
        spectra_ROI.coef_pol_x = coef_pol_x_ROI{roi} ;
        spectra_ROI.coef_pol_x(:,1) = spectra_ROI.coef_pol_x(:,1) + ...
            +delta_pos(2) ;
        spectra_ROI.coef_pol_dif = coef_pol_dif_ROI{roi} ;

        Model_spec = LinOpSpecProjLocLaw(sensor_ROI, ...
            spectra_ROI, pattern_model, false) * ...
            my_LinOpBroadcast([list_nb_lenslet(roi), nb_lambda], 1) ;

        %% Optimization
        C_spec = CostRobustPenalization( ...
            Model_spec, pic_ROI(:,:,roi), ...
            option_robust) ;

        % Regularization on the spectrum
        C_reg = CostL2([nb_lambda,1]) * ...
            LinOpGrad([nb_lambda,1], 1) ;


        % Global cost
        C = norm_spec*C_spec ;
        if mu_L2_grad>0
            C = C + norm_lamb*mu_L2_grad / ...
                (option_opti_spec.flag_s)^2*C_reg ;
        end

        % Optimization
        spec_roi(:,roi) = run_Opti(C, cons_min, cons_max, ...
            option_opti_spec, spec_roi(:,roi)) ;

        %% Display
        j_roi = 1+mod(roi-1, nb_tot_ROI^0.5) ;
        i_roi = 1+floor(mod((roi-1)/nb_tot_ROI^0.5, nb_tot_ROI^0.5)) ;

        j_roi = j_roi + (j_roi-1)*(2*rad_ROI_shift+1) + ...
            (1:(2*rad_ROI_shift+1)) ;
        i_roi = i_roi + (i_roi-1)*(2*rad_ROI_shift+1) + ...
            (1:(2*rad_ROI_shift+1)) ;

        pic_sim_ROI(i_roi,j_roi) = Model_spec*spec_roi(:,roi) ;
        pic_W_coef_ROI(i_roi,j_roi) = C_spec.computeW_(spec_roi(:,roi)) ;
    end
    save_fits(pic_sim_ROI, ['sim_spec_it_', num2str(cc)], ...
        save_path_shift_it) ;
    save_fits(pic_W_coef_ROI, ['W_coef_spec_it_', num2str(cc)], ...
        save_path_shift_it) ;
    save([save_path_shift_it, 'spec_it_', num2str(cc)], 'spec_roi', ...
        '-ascii', '-double') ;

    %% Estimation of the shift
    disp('   Estimating shift...') ;
    Model_shift = cell(nb_tot_ROI, 1) ;
    
    % Building the model
    for roi = 1:nb_tot_ROI
        Model_lenslet = cell(list_nb_lenslet(roi), 1) ;
        for lenslet_i = 1:list_nb_lenslet(roi)
            % Local law of the lenslet
            p_y = Pol_y * coef_pol_y_ROI{roi}(lenslet_i,:)' ;
            p_x = Pol_x * coef_pol_x_ROI{roi}(lenslet_i,:)' ;
            p_dif = Pol_dif * coef_pol_dif_ROI{roi}(lenslet_i,:)' ;
            
            % Elementary pattern
            Model_lambda = cell(nb_lambda, 1) ;
            for lambda = 1:nb_lambda
                Model_lambda{lambda} = get_pattern_operator( ...
                    sensor_y-p_y(lambda), sensor_x-p_x(lambda), [], [], ...
                        p_dif(lambda), spec_roi(lambda,roi), ...
                        0, pattern_model) ;
            end
            Model_lenslet{lenslet_i} = MapSummation(Model_lambda, ...
                ones(nb_lambda,1)) ;
        end
        Model_shift{roi} = MapSummation(Model_lenslet, ...
            ones(list_nb_lenslet(roi),1)) ;
    end

    % Building cost
    C_shift = cell(nb_tot_ROI, 1) ;
    for roi = 1:nb_tot_ROI
        % Option of the robust penalization
        option_robust = [] ;
        option_robust.method = option_opti_shift.RL2_method ;
        option_robust.noise_model = option_opti_shift.noise_model ;
        option_robust.var_0 = var_0 ;
        option_robust.eta = eta ;
        option_robust.flag_s = IFS_W_ROI(:,:,roi) .* ...
            option_opti_shift.flag_s ;
        
        C_shift{roi} = CostRobustPenalization( ...
            Model_shift{roi}, pic_ROI(:,:,roi), ...
            option_robust) ;
    end
    
    C = CostSummation(C_shift, ones(nb_tot_ROI,1)) ;
    
    % Optimization
    delta_pos = run_Opti(C, [], [], option_opti_shift, delta_pos) ;
    
    % Display
    pic_sim_ROI = zeros(nb_tot_ROI^0.5*(2*rad_ROI_shift+1) + ...
        nb_tot_ROI^0.5 +1 ) ;
    pic_W_coef_ROI = zeros(nb_tot_ROI^0.5*(2*rad_ROI_shift+1) + ...
        nb_tot_ROI^0.5 +1 ) ;
    for roi = 1:nb_tot_ROI
        j_roi = 1+mod(roi-1, nb_tot_ROI^0.5) ;
        i_roi = 1+floor(mod((roi-1)/nb_tot_ROI^0.5, nb_tot_ROI^0.5)) ;

        j_roi = j_roi + (j_roi-1)*(2*rad_ROI_shift+1) + ...
            (1:(2*rad_ROI_shift+1)) ;
        i_roi = i_roi + (i_roi-1)*(2*rad_ROI_shift+1) + ...
            (1:(2*rad_ROI_shift+1)) ;

        pic_sim_ROI(i_roi,j_roi) = Model_shift{roi}*delta_pos ;
        pic_W_coef_ROI(i_roi,j_roi) = ...
            C_shift{roi}.computeW_(delta_pos) ;
    end
    save_fits(pic_sim_ROI, ['sim_delta_it_', num2str(cc)], ...
        save_path_shift_it) ;
    save_fits(pic_W_coef_ROI, ['W_coef_delta_it_', num2str(cc)], ...
        save_path_shift_it) ;
    save([save_path_shift_it, 'delta_pos_it_', num2str(cc)], ...
        'delta_pos', '-ascii', '-double') ;
end

%% Saving
save([save_path_shift, 'delta_pos'], 'delta_pos', '-ascii', '-double') ;
close all ;