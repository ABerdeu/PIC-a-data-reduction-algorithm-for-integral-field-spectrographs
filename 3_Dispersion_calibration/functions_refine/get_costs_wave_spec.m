%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get the cost functions on the wavelength and spectral
% calibration files
% 
% Created: 03/21/2019 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #[coef_pol_y, coef_pol_x]# coefficients of the position polynomial laws
%
% #coef_pol_dif# coefficients of the diffraction polynomial law
%
% #list_lambda_cal# list of the wavelengths in the wavelength calibration
%
% #list_lambda# list of the wavelengths in the model
% 
% #trans# transmission of the lenslet
%
% #lamp_spec# spectrum of the lamp
%
% #pattern_model# pattern model to project
%
% #rad_ROI# radius of the region of interest to extract
%
% #rad_ROI_model# radius of the region of interest to simulate
%
% #pix# pixel characteristics of the sensor
%
% #[Calib_wave, Calib_spec]# calibration files
%
% #map_BP# Flag on the bad pixels
%
% #[RL2_method, noise_model, var_0, eta, flag_s]# Parameters for the robust
% penalization
%
%%%%%%%%
% Ouput
%%%%%%%%
% #[C_wave, C_spec]# the two cost functions
%
% #[Pol_y, Pol_x, Pol_dif]# the polynomials to simulate the spectral
% projection model
%
% #[crop_y, crop_x]# the extracted positions
%
% #Calib_spec# the extracted spectral calibration region
%
% #map_W_ij# the extracted weights
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C_wave, C_spec, Pol_y, Pol_x, Pol_dif, crop_y, crop_x, ...
    Calib_spec_ij, map_W_ij] = ...
    get_costs_wave_spec(coef_pol_y, coef_pol_x, coef_pol_dif, ...
    list_lambda_cal, list_lambda, trans, lamp_spec, pattern_model, ...
    pix, rad_ROI, rad_ROI_model, ...
    Calib_wave, Calib_spec, map_BP, ...
    RL2_method, noise_model, var_0, eta, flag_s)

    %% Degree of the polynomial expression
    deg_pol_y = length(coef_pol_y)-1 ;
    deg_pol_x = length(coef_pol_x)-1 ;
    deg_pol_dif = length(coef_pol_dif)-1 ;
    nb_lambda_cal = length(list_lambda_cal) ;
    nb_lambda = length(list_lambda) ;
    
    %% Extracting the region of interest
    [list_i, list_j] = ...
        get_ROI(list_lambda, coef_pol_y', ...
        coef_pol_x', rad_ROI, pix) ;

    [crop_y, crop_x] = get_extracted_pos(list_i, list_j, pix) ;

    %% Simulation of the residues
    spectra = [] ;
    spectra.list_lambda = list_lambda ;
    spectra.PSF_nb_x = 2*rad_ROI_model+1 ;
    spectra.PSF_nb_y = 2*rad_ROI_model+1 ;
    spectra.coef_pol_y = coef_pol_y ;
    spectra.coef_pol_x = coef_pol_x ;
    spectra.trans = trans ;
    spectra.coef_pol_dif = coef_pol_dif ;

    pix_spec = [] ;
    pix_spec.dy = 1 ;
    pix_spec.dx = 1 ;
    pix_spec.nb_y = length(list_i) ;
    pix_spec.nb_x = length(list_j) ;

    % Shift to the local frame
    spectra.coef_pol_y(1) = spectra.coef_pol_y(1) - ...
        crop_y(floor(pix_spec.nb_y/2+1)) ;
    spectra.coef_pol_x(1) = spectra.coef_pol_x(1) - ...
        crop_x(floor(pix_spec.nb_x/2+1)) ;

    % Local model
    Spec_proj_spec = LinOpSpecProjLocLaw(pix_spec, spectra, ...
        pattern_model, false) ;

    % Simulation
    Calib_spec_ij = extract_ij(Calib_spec, ...
        {list_i, list_j}, [pix.nb_y, pix.nb_x]) + ...
        Spec_proj_spec*lamp_spec' ;


    %% Optimizing the position laws

    % Model for the wavelength calibration

    % Estimation of the polynomial y-law
    Pol_y_cal_par = zeros(nb_lambda_cal, deg_pol_y+1) ;
    for d = 0:deg_pol_y
        Pol_y_cal_par(:,d+1) = list_lambda_cal.^d ;
    end
    Pol_y = zeros(nb_lambda, deg_pol_y+1) ;
    for d = 0:deg_pol_y
        Pol_y(:,d+1) = list_lambda.^d ;
    end

    % Estimation of the polynomial x-law
    Pol_x_cal_par = zeros(nb_lambda_cal, deg_pol_x+1) ;
    for d = 0:deg_pol_x
        Pol_x_cal_par(:,d+1) = list_lambda_cal.^d ;
    end
    Pol_x = zeros(nb_lambda, deg_pol_x+1) ;
    for d = 0:deg_pol_x
        Pol_x(:,d+1) = list_lambda.^d ;
    end

    % Estimation of the polynomial diffraction law
    Pol_dif_cal_par = zeros(nb_lambda_cal, deg_pol_dif+1) ;
    for d = 0:deg_pol_dif
        Pol_dif_cal_par(:,d+1) = list_lambda_cal.^d ;
    end
    Pol_dif = zeros(nb_lambda, deg_pol_dif+1) ;
    for d = 0:deg_pol_dif
        Pol_dif(:,d+1) = list_lambda.^d ;
    end

    % Local law of the lenslet
    p_dif = Pol_dif_cal_par * coef_pol_dif' ;

    % Elementary patterns for the spectral calibration
    Model_lambda = cell(nb_lambda_cal, 1) ;
    for lambda = 1:nb_lambda_cal
        sel_lambda = zeros(1, nb_lambda_cal) ;
        sel_lambda(lambda) = 1 ;
        Model_lambda{lambda} = get_pattern_operator( ...
            crop_y', crop_x', [], [], ...
                p_dif(lambda), ...
                [], ...
                0, pattern_model) * LinOpMatrix( ...
                [Pol_y_cal_par(lambda,:), ...
                    zeros(1, deg_pol_x+1 + nb_lambda_cal) ;...
                zeros(1, deg_pol_y+1), Pol_x_cal_par(lambda,:), ...
                    zeros(1, nb_lambda_cal) ; ...
                zeros(1, deg_pol_y+1 + deg_pol_x+1), sel_lambda]) ;
    end
    Model_wave = MapSummation(Model_lambda, ...
        ones(nb_lambda_cal,1)) ;

    % Model for the spectrum calibration

    % Selecting the spectrum calibration
    Sel_spec = true(deg_pol_x+1 + deg_pol_y+1 + nb_lambda_cal, 1) ;
    Sel_spec(deg_pol_x+1 + deg_pol_y+1 + (1:nb_lambda_cal)) = false ;
    Sel_spec = LinOpSelector(Sel_spec) ;

    % Local law of the lenslet
    p_dif = Pol_dif * coef_pol_dif' ;

    % Elementary patterns for the spectral calibration
    Model_lambda = cell(nb_lambda, 1) ;
    for lambda = 1:nb_lambda
        Model_lambda{lambda} = get_pattern_operator( ...
            crop_y', crop_x', [], [], ...
                p_dif(lambda), ...
                trans*lamp_spec(lambda), ...
                0, pattern_model) * LinOpMatrix( ...
                [Pol_y(lambda,:), zeros(1, deg_pol_x+1) ;...
                zeros(1, deg_pol_y+1), Pol_x(lambda,:)]) ;
    end
    Model_spec = MapSummation(Model_lambda, ...
        ones(nb_lambda,1))*Sel_spec ;

    % Weigth on the defective pixels
    map_W_ij = double(extract_ij(map_BP, {list_i, list_j}, ...
        [pix.nb_y, pix.nb_x])) ;
    map_W_ij(map_W_ij==0) = Inf ;

    % Option of the robust penalization and cost
    option_robust = [] ;
    option_robust.method = RL2_method ;
    option_robust.noise_model = noise_model ;
    option_robust.var_0 = var_0 ;
    option_robust.eta = eta ;
    option_robust.flag_s = map_W_ij .* flag_s ;

    C_wave = CostRobustPenalization( ...
        Model_wave, extract_ij(Calib_wave, ...
        {list_i, list_j}, ...
        [pix.nb_y, pix.nb_x]), option_robust) ;

    C_spec = CostRobustPenalization( ...
        Model_spec, Calib_spec_ij, option_robust) ;
end

