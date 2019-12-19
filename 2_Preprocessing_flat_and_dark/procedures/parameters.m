%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to set the parameters of the simulation
% 
% Created: 01/26/2019 (mm/dd/yyyy)
% Modified: 02/25/2019 (mm/dd/yyyy) Analyzing all dark current
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Data
% Path to the renamed raw data downloaded from the data center
data.path = 'D:\Anthony\Documents\2019_PostDoc Interlude\Data\2018_04_13 SAM_Faustine\Files_renamed/' ;

% Suffix of the full field sensor flats with the different lamps
data.flat_name = 'CALIB_FLAT_LAMP' ;

% Suffix of the dark acquisition with the shutted closed
data.dark_name = 'CALIB_DARK_OPEN' ;

% Suffix of the dark acquisition with the shutted opened
data.dark_BG_name = 'CALIB_DARK_BACKGROUND_OPEN' ;

% Suffix of the spectral flat acquisition (broadband source through the
% lenslet array)
data.spectral_flat_name = 'SPEC_FLAT_LAMP_PRIS_YJH' ;

% Suffix of the background sky acquisition
data.sky_BG_name = 'CALIB_DARK_BACKGROUND_OPEN' ; % Using the instrumenal
    % bakcground
% data.sky_BG_name = 'SCIENCE_SKY' ; % Using the acquisition on the sky

% Option of the spectral flat threholding
rad_dil = 7 ; % Radius to perform a morphological dilation
rad_med = 3 ; % Side of the pattern to perform a median filter 
threshold_data = 0.25 ; % Thershold to compute the binary mask to suppress
    % the background

tol_MAD = 15 ; % Tolerance in the dark for bad pixels indentification
size_med_filter = 5 ; % Size of the median filter to indentify
    % the bad pixels
    
nb_BP_edges = 4 ; % Number of not useful pixels on the edges of the sensor.
    
%% Plot parameters
fig_pos = [500, 500] ;
fig_size = [900, 600] ;
LineWidth = 3 ;
MarkerSize = 8 ;
FontSize = 25 ;
FontSize_axis = 15 ;

%% Saving parameters
save_path = './results/' ;
save_lin_fits = false ; % Flag to save the linear fits
