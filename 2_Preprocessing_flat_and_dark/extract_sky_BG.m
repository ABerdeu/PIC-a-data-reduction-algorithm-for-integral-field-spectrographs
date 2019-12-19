%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This macro is a procedure to extract the sky background from the remanent
% signal
%
% Created: 03/25/2019 (mm/dd/yyyy)
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

%% Set parameters
% Path the the results of the preprocessing of a given dataset
path_cal = './results/path_to_the_preprocessing_of_the_data_set/' ;

% Name of the last frame in data.path whose contribution must be removed
% from [data.path, data.sky_BG_name]
data.last_frame = 'name_of_the_last_frame.fits' ;

% How to deal the last frame sequence ?
%   'MEAN' -> mean of the sequence
%   'LAST' -> last frame of the sequence
flag_frame = 'MEAN' ;

% Number of loops to estimate the parameters
nb_cal = 100 ;


%% Building path
addpath(genpath(path_cal)) ;
copyfile('./extract_sky_BG.m', [path_cal, 'code/']) ;


%% Running parameters
run parameters

%% Opening files
% bad pixel map
BP_mask = fitsread([path_cal, 'IFS_BP_mask.fits']) ;
dark_BG = fitsread([path_cal, 'IFS_dark_BG.fits']) ;
BP_mask = BP_mask>0 ;

% last frame
name_f = [data.path, data.last_frame] ;
last_frame = fitsread(name_f) ;
exp_time_last_frame = get_fits_exposure(name_f) ;

switch flag_frame
    case 'MEAN'
        last_frame = mean(last_frame,3) ;
    case 'LAST'
        last_frame = last_frame(:,:,size(last_frame, 3)) ;
    otherwise
        error([flag_frame, ' is a unknown flag...']) ;
end

% Opening sky backgrounds
disp('Opening sky backgrounds...') ;
[nb_sky_BG, list_exp_time_sky_BG, list_avg_sky_BG, ...
    list_var_sky_BG] = ...
    open_list_files(data.path, data.sky_BG_name) ;
disp('Opening sky backgrounds: done!') ;

% Temporal frames
last_frame = last_frame/exp_time_last_frame ;
for t = 1:nb_sky_BG
    list_avg_sky_BG(:,:,t) = ...
        list_avg_sky_BG(:,:,t)/list_exp_time_sky_BG(t) ;
end
save_fits(last_frame, 'IFS_last_frame', path_cal) ;


%% Building model matrix
nb_pix = sum(BP_mask(:)) ;
vect_sky_BG = zeros(nb_pix, nb_sky_BG) ;
vect_data = last_frame(BP_mask)-dark_BG(BP_mask) ;

for f = 1:nb_sky_BG
    % Building the data vector
    sky_BG_f = list_avg_sky_BG(:,:,f) ;
    vect_sky_BG(:,f) = sky_BG_f(BP_mask)-dark_BG(BP_mask) ;
end

%% Optimization of the background and the remanent
sky_BG_b = mean(vect_sky_BG, 2) ;
sky_BG_a = vect_data - sky_BG_b ;
tau = ones(nb_sky_BG, 1) ;

[per_aux, lastprint] = display_percentage('init', ...
    '    Sky background identification') ;
delta_per = 5 ;
% Loop on the calibration
for cc = 1:nb_cal
    % Display the percentage
    [per_aux, lastprint] = display_percentage('iter', ...
        {cc/nb_cal*100, per_aux, delta_per, lastprint}) ;
    
    % Update of the remanent proportion in the sky background
    tau = ones(nb_sky_BG, 1) ;
    for f = 1:nb_sky_BG
        tau(f) = sum((vect_sky_BG(:,f) - sky_BG_b).*sky_BG_a) / ...
            sum(sky_BG_a.^2) ;
        
        % The ratio mus be decreasing
        tau(f) = min(tau) ;
    end
    tau = max(tau, 0) ;

    figure(1) ;
    hold on
    plot(tau) ;
    hold off
    pause(0.1) ;
    
    % Update of the remanent signal
    sky_BG_a = vect_data-sky_BG_b ;
    
    % Update of the sky background
    sky_BG_b = median(vect_sky_BG-sky_BG_a.*tau', 2) ;
    sky_BG_b = max(sky_BG_b,0) ;
    
    
    % Reshaping
    IFS_sky_BG_a = zeros(size(last_frame)) ;
    IFS_sky_BG_b = zeros(size(last_frame)) ;
    IFS_sky_BG_a(BP_mask) = sky_BG_a ;
    IFS_sky_BG_b(BP_mask) = sky_BG_b ;
    save_fits(IFS_sky_BG_a, ...
        ['IFS_sky_BG_a_', num2str(cc)], [path_cal, 'temp/IFS_sky_BG_a/']) ;
    save_fits(IFS_sky_BG_b, ...
        ['IFS_sky_BG_b_', num2str(cc)], [path_cal, 'temp/IFS_sky_BG_b/']) ;
end
% Cleaning percentage
display_percentage('exit', lastprint) ;

% Reshaping
IFS_sky_BG_a = zeros(size(last_frame)) ;
IFS_sky_BG_b = zeros(size(last_frame)) ;
IFS_sky_BG_a(BP_mask) = sky_BG_a ;
IFS_sky_BG_b(BP_mask) = sky_BG_b + dark_BG(BP_mask) ;
save_fits(IFS_sky_BG_a, 'IFS_sky_BG_a', path_cal) ;
save_fits(IFS_sky_BG_b, 'IFS_sky_BG_b', path_cal) ;