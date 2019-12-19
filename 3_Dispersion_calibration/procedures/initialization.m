%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to initialize the excution
% 
% Created: 10/10/2018 (mm/dd/yyyy)
% Modified: 02/06/2019 (mm/dd/yyyy) Modification to account for the
% background substraction
% Modified: 03/05/2019 (mm/dd/yyyy) Joint estimation with the spectral flat
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Opening sensor calibration
% Mask on the bad pixels
IFS_BP = fitsread([data.sensor_calib, 'IFS_BP_mask.fits']) ;

% Sensor flat
IFS_sensor_flat = fitsread([data.sensor_calib, 'IFS_sensor_flat.fits']) ;

% Temporal dark current
IFS_dark_current = fitsread([data.sensor_calib, 'IFS_dark_BG.fits']) ;

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

% Parameters of the Poisson noise
option_opti_amp_pos.eta = eta ;
option_opti_amp_pos.var_0 = var_0 ;
option_opti_par.eta = eta ;
option_opti_par.var_0 = var_0 ;
option_opti_spec.eta = eta ;
option_opti_spec.var_0 = var_0 ;
option_opti_spec_init.eta = eta ;
option_opti_spec_init.var_0 = var_0 ;


%% Opening data and correction by the sensor calibration
% Wave calibration file: correction by the dark current and the flat
IFS_calib_wave = fitsread(data.calib_wave) ;
IFS_calib_wave = mean(IFS_calib_wave,3) ; % Stacking the different frames
% % save_fits(IFS_calib_wave, 'IFS_calib_wave_not_corrected', save_path) ;
IFS_dark_wave = get_fits_exposure(data.calib_wave)*IFS_dark_current ./ ...
    IFS_sensor_flat ;
IFS_calib_wave = IFS_calib_wave./IFS_sensor_flat - IFS_dark_wave ;
IFS_BP(isnan(IFS_calib_wave)) = 0 ;
IFS_calib_wave(isnan(IFS_calib_wave)) = 0 ;
IFS_BP(IFS_calib_wave<=0) = 0 ;
IFS_calib_wave(IFS_calib_wave<=0) = 1 ;

save_fits(IFS_calib_wave, 'IFS_calib_wave_corrected', save_path) ;
save_fits(uint8(IFS_BP), 'IFS_BP_corrected', save_path) ;

% Wave calibration file: correction by the dark current and the flat
IFS_calib_spec = fitsread(data.calib_spec) ;
IFS_calib_spec = mean(IFS_calib_spec,3) ; % Stacking the different frames
% % save_fits(IFS_calib_spec, 'IFS_calib_spec_not_corrected', save_path) ;
IFS_dark_spec = get_fits_exposure(data.calib_spec)*IFS_dark_current ./ ...
    IFS_sensor_flat ;
IFS_calib_spec = IFS_calib_spec./IFS_sensor_flat - IFS_dark_spec ; 
IFS_BP(isnan(IFS_calib_spec)) = 0 ;
IFS_calib_spec(isnan(IFS_calib_spec)) = 0 ;
IFS_BP(IFS_calib_spec<=0) = 0 ;
IFS_calib_spec(IFS_calib_spec<=0) = 1 ;

save_fits(IFS_calib_spec, 'IFS_calib_spec_corrected', save_path) ;

% Pixel characteristic
[pix_IFS.nb_y, pix_IFS.nb_x, pix_IFS.nb_f] = size(IFS_calib_spec) ;
pix_IFS.dx = 1 ;
pix_IFS.dy = 1 ;

% Weigth on the defective pixels
IFS_W = IFS_BP ;
IFS_W(IFS_W==0) = Inf ;


%% Display dynamics
val_min_wave = 0 ;
val_max_wave = IFS_calib_wave((IFS_BP>0)&(IFS_mask>0)) ;
val_max_wave = 10*median(val_max_wave(:)) ;
val_min_spec = 0 ;
val_max_spec = IFS_calib_spec((IFS_BP>0)&(IFS_mask>0)) ;
val_max_spec = 3*median(val_max_spec(:)) ;

%% Spectra parameters
% Number of hexagons and wavelengths
if exist('range_lambda', 'var') && exist('delta_lambda', 'var')
    list_lambda = ...
        (min(range_lambda):delta_lambda:max(range_lambda))' ;
elseif exist('list_lambda_file', 'var')
    list_lambda = fitsread(list_lambda_file) ;
    list_lambda = fit_pol(deg_pol_list_lambda, ...
        ind_list_lambda', list_lambda', ind_list_lambda_new') ;
else
    error('There is no information on the lambda discretization...') ;
end
list_lambda = list_lambda(:) ; % Vector shape
save([save_path, 'list_lambda.txt'], 'list_lambda', '-ascii', '-double') ;
save([save_path, 'list_lambda_cal.txt'], ...
    'list_lambda_cal', '-ascii', '-double') ;
nb_lambda = length(list_lambda) ;

% Changing the origin the first calibration wavelength (insure to control
% the divergence of the parameters
lambda_0 = mean(list_lambda_cal(:)) ;
save([save_path, 'lambda_0.txt'], 'lambda_0', '-ascii', '-double') ;
list_lambda = list_lambda-lambda_0 ;
list_lambda_cal = list_lambda_cal- lambda_0 ;

%% Initialization
nb_lambda_cal = length(list_lambda_cal) ; % Number of wavelength
list_color = jet(nb_lambda_cal) ;

fig_1 = figure(1) ;
subplot(1,2,1) ;
imshow(IFS_calib_wave, [val_min_wave, val_max_wave]) ;
title('Click on an area to zoom in...') ;
subplot(1,2,2) ;
imshow(IFS_calib_spec, [val_min_spec, val_max_spec]) ;
title('Click on an area to zoom in...') ;

% User pointing
p_center = zeros(1,2) ;
[p_center(2), p_center(1)] = my_ginput(1, 0.5*[1,1,1]) ;
p_center = round(p_center) ;

% Insure the extracted frame is in picture
nb_center = round((per_cen/100).^0.5*[pix_IFS.nb_y, pix_IFS.nb_x]) ;
nb_center = 2*round(nb_center/2) ;
p_center(2) = max(nb_center(2)/2+1, min(pix_IFS.nb_y-nb_center(2)/2, ...
    p_center(2))) ;
p_center(1) = max(nb_center(1)/2+1, min(pix_IFS.nb_x-nb_center(1)/2, ...
    p_center(1))) ;

% Index to extract
y_center = get_Fourier_vector(nb_center(1),1)+p_center(1) ;
x_center = get_Fourier_vector(nb_center(2),1)+p_center(2) ;

p_center(2) = p_center(2) - floor(pix_IFS.nb_y/2+1) ;
p_center(1) = p_center(1) - floor(pix_IFS.nb_x/2+1) ;

IFS_zoom_wave = IFS_calib_wave(y_center, x_center) ;
IFS_zoom_spec = IFS_calib_spec(y_center, x_center) ;
pix_zoom.nb_y = nb_center(1) ;
pix_zoom.nb_x = nb_center(2) ;
pix_zoom.nb_f = 1 ;
pix_zoom.dx = 1 ;
pix_zoom.dy = 1 ;
close(fig_1) ;

%% Highlighting bad pixels
sel_BP = IFS_BP(y_center, x_center)==0 ;
IFS_zoom_wave_GB = IFS_zoom_wave ;
IFS_zoom_wave_GB(sel_BP) = 0.5*(val_max_wave-val_min_wave)+val_min_wave  ;
IFS_zoom_wave(sel_BP) = 1*(val_max_wave-val_min_wave)+val_min_wave  ;
IFS_zoom_wave = cat(3, IFS_zoom_wave, IFS_zoom_wave_GB, IFS_zoom_wave_GB) ;
IFS_zoom_wave = max(0, ...
    min(val_max_wave-val_min_wave, IFS_zoom_wave-val_min_wave)) / ...
    (val_max_wave-val_min_wave) ;

IFS_zoom_spec_GB = IFS_zoom_spec ;
IFS_zoom_spec_GB(sel_BP) = 0.5*(val_max_spec-val_min_spec)+val_min_spec ;
IFS_zoom_spec(sel_BP) = 1*(val_max_spec-val_min_spec)+val_min_spec  ;
IFS_zoom_spec = cat(3, IFS_zoom_spec, IFS_zoom_spec_GB, IFS_zoom_spec_GB) ;
IFS_zoom_spec = max(0, ...
    min(val_max_spec-val_min_spec, IFS_zoom_spec-val_min_spec)) / ...
    (val_max_spec-val_min_spec) ;

clear('sel_BP', 'IFS_zoom_spec_GB', 'IFS_zoom_wave_GB') ;
