%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to initialize the reconstruction
% 
% Created: 06/01/2018 (mm/dd/yyyy)
% Modified: 08/08/2018 (mm/dd/yyyy) Taking into account the transmission
% given by the flat. Reweighted L2 norm
% Modified: 10/16/2018 (mm/dd/yyyy) Taking into account the full field flat
% and the Poisson noise in the robust penalization
% Modified: 11/20/2018 (mm/dd/yyyy) Forcing edges of the spectrum to be
% null and new regularizations
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Initialization...') ;

disp('   Loading calibration and data...') ;

%% Loading sensor calibration
% Mask on the bad pixels
IFS_BP = fitsread([path.calib_sensor, 'IFS_BP_mask.fits']) ;

% Sensor flat
IFS_sensor_flat = fitsread([path.calib_sensor, 'IFS_sensor_flat.fits']) ;

% Temporal dark current
IFS_BG = fitsread([path.calib_sensor, 'IFS_', flag_BG, '.fits']) ;

% Mask on the scientific area
IFS_mask = fitsread([path.calib_sensor, 'IFS_mask.fits']) ;

% Variance law for the Poisson noise
if flag_MAD
    eta = load([path.calib_sensor, 'eta_MAD.txt'], '-ascii') ;
    var_0 = load([path.calib_sensor, 'sig_MAD.txt'], '-ascii') ;
else
    eta = load([path.calib_sensor, 'eta.txt'], '-ascii') ;
    var_0 = load([path.calib_sensor, 'sig.txt'], '-ascii') ;
end

%% Preprocessing the data
IFS_data = fitsread([data.path data.name, '.fits']) ;
IFS_data = mean(IFS_data,3) ; % Stacking the different frames
IFS_BG = exp_time*IFS_BG ./ IFS_sensor_flat ;
IFS_data = IFS_data./IFS_sensor_flat ;
IFS_BP(isnan(IFS_data)) = 0 ;
IFS_data(isnan(IFS_data)) = 0 ;
IFS_BP(IFS_data<=0) = 0 ;
IFS_data(IFS_data<=0) = 1 ;

save_fits(IFS_data, 'IFS_data_corrected', save_path) ;
save_fits(IFS_data-IFS_BG, 'IFS_data_corrected_BG', save_path) ;
save_fits(IFS_BG, 'IFS_BG', save_path) ;
save_fits(IFS_BP, 'IFS_BP', save_path) ;

% Parameters of the IFS
[pix_IFS.nb_y, pix_IFS.nb_x] = size(IFS_data) ;
pix_IFS.dx = 1 ;
pix_IFS.dy = 1 ;


% Wright on the defective pixels
IFS_BP_W = ones([pix_IFS.nb_y, pix_IFS.nb_x]) ;
IFS_BP_W(IFS_BP==0) = Inf ;
save_fits(uint8(1./IFS_BP_W), 'IFS_BP/IFS_BP_it_0', save_path) ;

disp(   'Loading done!') ;


%% Model
disp('   Loading model...') ;
global Conv_hexa ;  load([data.model_path, 'Conv_hexa.mat']) ;
global Interp ;     load([data.model_path, 'Interp.mat']) ;
global Pad ;        load([data.model_path, 'Pad.mat']) ;
global Spec_proj ;	load([data.model_path, 'Spec_proj.mat']) ;

disp(   'Loading done!') ;

%% Parameters
% Parameters of the hypercube
pix_hypcube.nb_y = Pad.sizein(1) ;
pix_hypcube.nb_x = Pad.sizein(2) ;
pix_hypcube.nb_l = Pad.sizein(3) ;

nb_lenslet = Spec_proj.sizein(1) ;
nb_lambda = Spec_proj.sizein(2) ;

if isfield(data, 'pos_center')
    if ~isempty(data.pos_center)
        warning(['data.pos_center is existing... ', ...
            'Redefinition of the field center...']) ;
    end
    pos_center = load(data.pos_center, '-ascii') ;
    if mod(pix_hypcube.nb_y,2)==0 % The zero position for SPHERE is between
        % pixels
        pos_center(1) = pos_center(1)+0.5 ;
    end
    if mod(pix_hypcube.nb_x,2)==0 % The zero position for SPHERE is between
        % pixels
        pos_center(2) = pos_center(2)+0.5 ;
    end
    list_y = load([data.model_path, 'coef_model/coef_pol_y.txt'], ...
        '-ascii') ;
    list_x = load([data.model_path, 'coef_model/coef_pol_x.txt'], ...
        '-ascii') ;
    list_y = list_y(:,1) ;
    list_x = list_x(:,1) ;
    
    Interp = LinOpInterp({get_Fourier_vector(pix_hypcube.nb_y, 1)', ...
        get_Fourier_vector(pix_hypcube.nb_x, 1)'}, ...
        [list_y/pix_hypcube.dy-pos_center(1), ...
        list_x/pix_hypcube.dx-pos_center(2)], [], ...
        [pix_hypcube.nb_y, pix_hypcube.nb_x, nb_lambda], [1,2]) ;
end

%% Forward model
Forward_Mod = Spec_proj * Interp * Pad' * Conv_hexa * Pad ;

%% Declaration of the regularization costs
% Cost_TV
Cost_TV = CostHyperBolic( ...
    [pix_hypcube.nb_y, pix_hypcube.nb_x, pix_hypcube.nb_l, 2], ...
    eps,  3:4) * LinOpGrad( ...
    [pix_hypcube.nb_y, pix_hypcube.nb_x, pix_hypcube.nb_l], 1:2) ;


% Cost_L2_grad_lambda
Cost_L2_grad_lambda = CostL2( ...
    [pix_hypcube.nb_y, pix_hypcube.nb_x, pix_hypcube.nb_l]) * ...
    LinOpGrad( ...
    [pix_hypcube.nb_y, pix_hypcube.nb_x, pix_hypcube.nb_l], 3, ...
    'zeros') ;

% Mixed norm on the wavelength and the positions
% Cost_lambda = mu_lambda/(nb_y*nb_x*nb_l) * ...
%     CostHyperBolic([nb_y, nb_x, nb_l], eps, 3) ;

%% Declaration of the constraints
hyp_cub_min = zeros(pix_hypcube.nb_y, pix_hypcube.nb_x, pix_hypcube.nb_l) ;
hyp_cub_max = Inf(pix_hypcube.nb_y, pix_hypcube.nb_x, pix_hypcube.nb_l) ;
hyp_cub_max(:, :, [1, pix_hypcube.nb_l]) = max_val_edge ;

%% Clearing unused operators
clear('Spec_proj',  'Interp', ...
    'Pad', 'Conv_hexa') ;

disp('Initialization done!') ;