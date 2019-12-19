%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to set the parameters of the simulation
% 
% Created: 05/31/2018 (mm/dd/yyyy)
% Modified: 08/08/2018 (mm/dd/yyyy) New functions
% Modified: 10/25/2018 (mm/dd/yyyy) Taking into account the spectral flat
% fitting of the parmeters
% Created: 11/27/2018 (mm/dd/yyyy) Adaptation to polynomial laws
% Created: 03/13/2019 (mm/dd/yyyy) Adaptation to local polynomial laws and
% autocalibration on the data
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Data
% Path to the calibration refinement
path_calib_spec = ['../3_Dispersion_calibration/results/', ...
    'path_to_the_calibration/refinement_pos_dif_par/coef/'] ;

% Calibration number
c_cal = 10 ;


%% Shift determination
% Need to perform an autocalibration on the science data?
flag_estimate_shift = true ;

% Exposure time found in the header of the fits file of the data
exp_time = 4 ;

% Path and name of the science data on which the model must be aligned
path_data = 'path_to_the_science_data.fits' ;

% Path to the preprocessing and the sensor calibration
path_calib_sensor = ['../2_Preprocessing_flat_and_dark/results/', ...
    'path_to_the_preprocessing/'] ;

% Flag to select the variance law (estimated via MAD or not)
flag_MAD = true ; 

% Number of calibration aleterning the spectrum fit and the shift
% refinement
nb_cal = 3 ;

% Parameters of the extracted regions of interest
nb_ROI = 3 ; % Number of extracted regions of interest along a side
theta_square = -10 ; % Orientation of the selected ROI
rad_ROI_shift = 30 ; % Radius of the region of interest to extract

% Option of the optimizer for the positions
option_opti_shift.RL2_method = 'Cauchy' ;
option_opti_shift.noise_model = 'none' ;
option_opti_shift.method = 'fminsearch' ;
option_opti_shift.verbose = false ;
option_opti_shift.maxiter = 100 ; % Maximal number of iterations for the
option_opti_shift.flag_s = 350 ;

% Option of the optimizer for the spectra indentification
option_opti_spec = option_opti_shift ;
option_opti_spec.method = 'VMLMB' ;

% Maximal value on the spectral edge
max_val_edge = 0*1.5e4 ;

% Regularization on the spectrum
mu_L2_grad = 8e-4 ;

%% Parameters of the cartesian grid on the lenslet
%   - pix.[nb_x, nb_y] -> The size of the domain
%   - pix.[nb_x_pad, nb_y_pad] -> The size of the padded domain (if 0, no
%   padding)
%   - pix.[dx, dy] -> The pixel side pitch

% Pixel size (must be compared to the sensor size)
pix.dx = 8 ;
pix.dy = 8 ;

% Number of pixels
pix.nb_x = 300 ;
pix.nb_y = 300 ;

% Padding (empty = no padding)
pix.nb_x_pad = 2*pix.nb_x ;
pix.nb_y_pad = 2*pix.nb_y ;


%% Parameters of the cartesian grid on the sensor
%   - sensor.[nb_x, nb_y] -> The size of the domain
%   - sensor.[dx, dy] -> The pixel side pitch
%   - sensor.[PSF_nb_x, PSF_nb_y] -> Spatial extension of the 
%   elementary PSF for a given lenslet and for a given wavelength (in
%   pixel)

% Pixel size
sensor.dx = 1 ;
sensor.dy = 1 ;

% Number of pixels
sensor.nb_x = 2048 ;
sensor.nb_y = 2048 ;

%% Parameters for the lenslet pattern
% Pattern model (axisymmetric)
pattern_model.flag_norm = true ; % Normalized?
pattern_model.oversampling = 1 ; % Scale to oversample the patternt to fit 
    % (must be an integer)
% Profile
%   Gaussian -> Gaussian pattern (sigma)
%   Moffat -> Moffat pattern (alpha, beta)
pattern_model.flag_profile = 'Gaussian' ;
% pattern_model.flag_profile = 'Moffat' ;

%% Parameters to simulate the lenslets optical properties
rad_ROI = 4 ; % Radius of the region of interest to simulate the PSF

% Fill factor of the lenslet (in terms of percentage of the side)
lenslet_fill_factor = 100 ; % (%)
    
% #flag_space#
%   SD -> compute the convolution kernel in the spatial domain
%   FD -> compte the convolution kernel in the Fourier domain
flag_space = 'SD' ;

% Numerical constant for the convolution operator definition
%   flag_space = SD -> factor of scaling to compute the apodization on the
%   edges of the hexagon
%   flag_space = FD -> numerical threshold to use Taylor expansion of the
%       formula (if 0, the numerical computation is supposed to be valid)
const = 500 ;


%% Saving parameters
save_path = './results/' ;

%% Plot parameters
LineWidth = 2 ;
MarkerSize = 10 * LineWidth ;
fig_pos = [75, 100] ;
fig_size = 900*[1, 1] ;
FontSize_axis = 15 ;
FontSize = 30 ;
