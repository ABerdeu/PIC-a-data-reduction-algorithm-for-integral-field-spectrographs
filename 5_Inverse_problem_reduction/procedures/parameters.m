%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to set the parameters of the simulation
% 
% Created: 06/01/2018 (mm/dd/yyyy)
% Modified: 08/08/2018 (mm/dd/yyyy) Taking into account the transmission
% given by the flat. Reweighted L2 norm
% Modified: 10/16/2018 (mm/dd/yyyy) Taking into account the full field flat
% and the Poisson noise in the robust penalization
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Data
% Path to the data and name of the file
data.path = '../1_Data/path_to_the_data_folder/' ;
data.name = 'file_name.fits' ;

% Path to the model
data.model_path = ['../4_PIC_operators_declaration/', ...
    'results/path_to_the_model/'] ;

% ¨Path to the sensor calibration and the files preprocessing
path.calib_sensor = ['../2_Preprocessing_flat_and_dark/results/', ...
    'path_to_the_preprocessing/'] ;

% Flag to select the variance law (estimated via MAD or not)
flag_MAD = true ;

% Flag to select the dark current
    % 'dark_BG' -> dark current with the shutter open but not on the sky
    % 'sky_BG'  -> dark current on the sky (careful with the remanence of
    % the sensor)
    % 'sky_BG_b'-> dark current on the sky (partially corrected for the 
    % remanence)
flag_BG = 'dark_BG' ;

% Exposure time found in the header of the fits file of the data
exp_time = 1.6507300 ;

% % Path to the position of the field center
% % (Comment if no specific center)
% data.pos_center = ['./results/', ...
%     'path_to_the_data_OBJECT_CENTER/', ...
%     'centroid_analysis/pos_center.txt'] ;
% pix_hypcube.dx = 8 ;
% pix_hypcube.dy = 8 ;

%% Optimizer parameters
% Number of residues estimation
nb_res_estimation = 6 ;


% Maximal value on the spectral edge
max_val_edge = 0*1.5e4 ;

% Hyperparameters
eps = 1e-16 ; % Epsilon for the computation of the norms
list_mu_TV = 5e-2*ones(nb_res_estimation, 1) ;
list_L2_grad_lambda = 0.5e-5*ones(nb_res_estimation, 1) ;

% Size of the median filter for the hypercube and the IFS
med_filter_IFS = [3, 3] ;
med_filter_cube = [5, 5, 3] ;


%% Optimizer parameters
% Method for the reweighted least square
%   ->  'Andrews'
%   ->  'Beaton-Tukey'
%   ->  'Cauchy'
%   ->  'Fair'
%   ->  'Huber'
%   ->  'Logistic'
%   ->  'Talwar-Hinich'
%   ->  'Welsch-Dennis'
option_opti.RL2_method = 'Cauchy' ;

% List of the scaling factor to scale the residues
list_scaling_factor = 25*ones(1,nb_res_estimation) ;

% Median filter on the reconstructed volume to erase the bad pixels
list_median_flag = false(1,nb_res_estimation) ;

% Is the Poisson noise dynamic in the cost function or not?
list_Poisson_flag = false(1,nb_res_estimation) ;

% Thresholding the weighting coefficients?
list_threshold_W_coef = true(1,nb_res_estimation) ;
list_threshold_W_coef(1) = false ;
th_W_coef = 0.5 ;

% Option of the optimizer
option_opti.method = 'VMLMB' ;
option_opti.verbose = true ;
option_opti.maxiter = 20 ; % Maximal number of iterations for the
    % optimizer
option_opti.memoizeOpts = true ; % Memorization of the costs computations
    % to save time in the iteration process

    
% Option of the VMLMB optimizer
option_opti.ItUpOut = 1 ;
option_opti.m = 3 ;  % number of memorized step in hessian approximation

%% Saving parameters
save_path = './results/' ;