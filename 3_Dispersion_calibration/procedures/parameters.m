%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to set the parameters of the simulation
% 
% Created: 04/27/2018 (mm/dd/yyyy)
% Modified: 07/24/2018 (mm/dd/yyyy) Adding automatic background suppression
% Modified: 08/31/2018 (mm/dd/yyyy) Modification to axisymmetric profiles
% Modified: 10/10/2018 (mm/dd/yyyy) Modification to identical profiles at
% each wavelength
% Modified: 02/06/2019 (mm/dd/yyyy) Modification to account for the
% background substraction
% Modified: 03/05/2019 (mm/dd/yyyy) Joint estimation with the spectral flat
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Data
% Wavelength calibration file
data.calib_wave = ['path_to_the_calibration_file_', ...
    'CALIB_WAVE_LAMP_PRIS_XXX.fits'] ;

% Broadband lamp spectral calibration file
data.calib_spec = ['path_to_the_calibration_file_', ...
    'CALIB_SPEC_FLAT_LAMP_PRIS_XXX.fits'] ;

% Path the the preprocessing of the dataset
data.sensor_calib = ['../2_Preprocessing_flat_and_dark/results/', ...
    'path_to_the_preprocessing/'] ;

% Flag to select the variance law (estimated via MAD or not)
flag_MAD = true ; 

flag_PRISM = 'YH' ; % 'YJ' / 'YH'

%% Region of interest to analyze the PSFs
rad_ROI = 2 ; % Radius of the region of interest to extract

% Lambda
switch flag_PRISM
    case 'YJ'
        %%%%%% Y_J mode -> range: 0.95 - 1.35 µm // calibration: 0.98772, 1.12371, 1.30937 %%%%%%
        list_lambda_cal = [0.98772, 1.12371, 1.30937]' ;
        
        % List of the model wavelengths
        list_lambda_file = './data/Lambda_info_Y_J.fits' ; % 1:39
        ind_list_lambda = 1:39 ;
        deg_pol_list_lambda = 3 ;
        ind_list_lambda_new = -3:1:45 ;
    case 'YH'
        %%%%%% Y_H mode -> range: 0.95 - 1.65 µm // calibration: 0.98772, 1.12371, 1.30937, 1.5451 %%%%%%
        % List of the calibration wavelengths
        list_lambda_cal = [0.98772, 1.12371, 1.30937, 1.5451]' ;
        
        % List of the model wavelengths
        list_lambda_file = './data/Lambda_info_H.fits' ; % 1:39
        ind_list_lambda = 1:39 ;
        deg_pol_list_lambda = 3 ;
        ind_list_lambda_new = -2:1:43 ;
        
    otherwise
        error([flag_PRISM, ': unknown flag for the prism (YJ/YH)']) ;
end

% Pattern model (axisymmetric)
pattern_model.flag_norm = true ; % Normalized?
pattern_model.flag_scale = false ; % Is the pattern scaled according to
    % lambda?
% pattern_model.scale = list_lambda_cal/list_lambda_cal(1) ; % Scaling rule
% pattern_model.scale = ones(size(list_lambda_cal)) ; % Scaling rule
pattern_model.oversampling = 1 ; % Scale to oversample the pattern to fit 
    % (must be an integer)
% Profile
%   Gaussian -> Gaussian pattern (sigma)
%   Hexagon -> Gaussian pattern (sigma)
%   Moffat -> Moffat pattern (alpha, beta)
pattern_model.flag_profile = 'Gaussian' ;
pattern_user_guess = 1 ;
% pattern_model.flag_profile = 'Hexagon' ;
% pattern_model.theta = 19 ;
% pattern_user_guess = 6 ;
% pattern_model.flag_profile = 'Moffat' ;
% pattern_user_guess = [2, 3] ; % Moffat

%% Calibration intialization
per_cen = 0.5 ; % Percentage of the field to display for the user manual
    % selection
% Number of calibration iteration to refine the first hexagon
% -> 1 iteration = gaussian fitting + centroid determination
nb_c = 10 ;
nb_spec_init = 19 ;

%% Robust penalization parameters
% Method for the reweighted least square
%   ->  'Andrews'
%   ->  'Beaton-Tukey'
%   ->  'Cauchy'
%   ->  'Fair'
%   ->  'Huber'
%   ->  'Logistic'
%   ->  'Talwar-Hinich'
%   ->  'Welsch-Dennis'
option_opti_amp_pos.RL2_method = 'Cauchy' ;
option_opti_amp_pos.noise_model = 'Poisson' ;

% Method to compute the residue scaling
%   ->  numeric -> constant value
%   ->  MAD     -> median absolute deviation of the region of interest
option_opti_amp_pos.flag_s = 25 ;


%% Polynomial fit
switch flag_PRISM
    case 'YJ'
        deg_pol_y = 2 ;
        deg_pol_x = 2 ;
        deg_pol_dif = 2 ;
        
    case 'YH'
        deg_pol_y = 3 ;
        deg_pol_x = 3 ;
        deg_pol_dif = 2 ;
end


%% Optimization parameters for the spectrum fit
% Calibration technique: which are the steps and the order of the
% operations realized for each haxagon?
% Examples:
%   {'pos_amp'} ; -> fitting both the position and the amplitudes of the
%   patterns
%   {'amp', 'pos'} ; -> fitting the amplitude and then the positions of the
%   patterns
%   'wave' -> fitting the wave calibration
%   'spec' -> fitting the spectrum only
%   'both' -> joint fit on the wave and the spectral clibration
list_cal_pat_init = {'amp', 'pos'} ; % Calibration technique for the
    % pattern fitting in the initialization
list_cal_spec_init = {'wave', 'spec', 'both', 'spec', 'both', 'spec', ...
    'both', 'spec', 'both'} ; % Must finish by 'both'
list_cal_spec_esti = {'spec', 'both', 'spec', 'both', 'spec', ...
    'both', 'spec', 'both'} ; % Must finish by 'both'

% Maximal value on the spectral edge
max_val_edge = 0*1.5e4 ;

% Weight on the wavelength calibration compared to the spectrum fit in the
% cost function
mu_wave = 0.5 ;

% Regularization on the spectrum
mu_L2_grad = 1.5e-3 ;

list_cal = {'pos_amp'} ;
tol_trans = 0.2 ; % Tolerance on the transmission to suppose that a
    % spectrum is correctly fitted




%% Optimization parameters
% Option of the optimizer for the position and amplitude in the
% initialization
option_opti_amp_pos.method = 'fminsearch' ;
option_opti_amp_pos.verbose = false ;
option_opti_amp_pos.maxiter = 250 ; % Maximal number of iterations for the
    % optimizer
option_opti_amp_pos.memoizeOpts = true ; % Memorization of the costs
    % computations to save time in the iteration process

% Option of the VMLMB optimizer (if used)
option_opti_amp_pos.ItUpOut = 1 ;
option_opti_amp_pos.m = 3 ;  % number of memorized step in hessian
    % approximation

% Option of the optimizer for the pattern parameters in the initialization
option_opti_par = option_opti_amp_pos ;
option_opti_par.maxiter = 250 ;
option_opti_par.verbose = false ;
option_opti_par.method = 'fminsearch' ;

% Option to fit the spectra
option_opti_spec_init = option_opti_par ;
option_opti_spec_init.noise_model = 'none' ;
option_opti_spec_init.flag_s = 350 ;
option_opti_spec_init.maxiter = 500 ;
option_opti_spec_init.verbose = false ;
option_opti_spec_init.method = 'VMLMB' ;
option_opti_spec = option_opti_spec_init ;
option_opti_spec.maxiter = 250 ;



%% Plotting option
LineWidth = 2 ;
    
%% Saving parameters
save_path = './results/' ;


