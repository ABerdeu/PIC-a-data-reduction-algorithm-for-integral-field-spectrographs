%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This macro is a procedure to automatically calibrate the lenslet array
% projection on the camera field
% /!\ If the calibration fits contains several frame, it is averaged.
%
% Created: 05/17/2018 (mm/dd/yyyy)
% Modified: 07/24/2018 (mm/dd/yyyy) Adding automatic background suppression
% Modified: 08/31/2018 (mm/dd/yyyy) Modification to axisymmetric profiles
% Modified: 10/10/2018 (mm/dd/yyyy) Modification to identical profiles at
% each wavelength
% Modified: 02/06/2019 (mm/dd/yyyy) Modification to account for the
% background substraction
% Modified: 03/05/2019 (mm/dd/yyyy) Joint estimation with the spectral flat
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
% Modified: 10/14/2019 (mm/dd/yyyy) On-line deposit
% Author: Anthony Berdeu (Centre de Recherche en Astrophysique de Lyon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% #spec# index of the spectrum
% #nb_spec# number of spectra
% #l# index of the wavelength (lambda)
% #nb_lambda_cal# number of wavelengths
% #coef_pol_x/coef_pol_y# list of the positions (h,l) in the field (even
%   outside the lenslet grid)
% #map_spec# map of the index of the spectra in the picture
%   -> -1: the pixel is out of the field of view
%   -> 0: the pixel is not associated with a known spectrum
%   -> i: index of the closest spectrum
%
clear all ;
close all ;
clc ;   

% Absolute path
abs_path = pwd ;
abs_path = [abs_path, '\'] ;

% Load functions and packages
restoredefaultpath
addpath(genpath('./functions')) ;
addpath('./procedures') ;

% Loading toolboxes
addpath(genpath('../Add-ons/')) ;


%% Set parameters
run parameters
save_path = [save_path, ...
    datestr(now,'yyyy_mm_dd_HH_MM_SS'),'/'] ;
save_path = make_dir(save_path) ;

% Saving the running codes
save_path_aux = make_dir([save_path, 'code/']) ;
copyfile('./executable.m', save_path_aux) ;
save_path_aux = make_dir([save_path, 'code/procedures/']) ;
copyfile('./procedures',save_path_aux) ;
save_path_aux = make_dir([save_path, 'code/functions/']) ;
copyfile('./functions',save_path_aux) ;


%% Initialization
run initialization


%% First spectrum analysis
first_spec = true ;
while first_spec
    run fit_first_spectrum
    first_spec = ~y_user ;
end

%% Robust estimation of the spectrum by fitting a local shared spectrum
run fit_lamp_spectrum
save([save_path, 'mu_wave_norm.txt'], 'mu_wave_norm', ...
    '-ascii', '-double') ;

%% Hexagonal grid parameters
side_hex = ...
    (coef_pol_x_init(2:7,1) - coef_pol_x_init(1,1)) + ...
    (coef_pol_y_init(2:7,1) - coef_pol_y_init(1,1))*1i ;
side_hex = mean(side_hex,2) ;
side_hex = mean( abs(side_hex).* ...
    exp(1i*mod(angle(side_hex)*180/pi,60)*pi/180)) ;
theta_hex = mod(180/pi*angle(side_hex)+30,60) ;
pattern_model.theta = theta_hex ;
side_hex = abs(side_hex) ;
save([save_path, '/theta_hex.txt'], 'theta_hex', '-ascii', '-double') ;
theta_hex = theta_hex-30 ;

%% Maximal number of hexagons and hexagon map
nb_rad = floor(side_hex/(2*rad_ROI)) ;
map_spec = zeros(pix_IFS.nb_y+2*round(side_hex), ...
    pix_IFS.nb_x+2*round(side_hex), ...
    'int16') ;
[pix_map.nb_y, pix_map.nb_x] = size(map_spec) ;
pix_map.dx = 1 ;
pix_map.dy = 1 ;

surf_hex = sum(IFS_mask(:)) ;


%% Initialization with the first spectrum
nb_h_max = floor(pix_map.nb_y*pix_map.nb_x/(nb_rad*rad_ROI*2+1))+1 ;
coef_pol_y = zeros(nb_h_max, deg_pol_y+1) ;
coef_pol_x = zeros(nb_h_max, deg_pol_x+1) ;
coef_pol_dif = zeros(nb_h_max, deg_pol_dif+1) ;
list_amp = zeros(nb_h_max, nb_lambda_cal) ;
trans = zeros(nb_h_max, 1) ;
IFS_sim_wave = zeros(pix_IFS.nb_y, pix_IFS.nb_x) ;
IFS_sim_spec = zeros(pix_IFS.nb_y, pix_IFS.nb_x) ;

% Refining positions on the wavelength and spectral calibration
[coef_pol_y(1,:), coef_pol_x(1,:), coef_pol_dif(1,:), ...
    trans(1), list_amp(1,:), list_i, list_j, Sim_wave, Sim_spec] = ...
    fit_calib_spectrum( ...
    coef_pol_y_init(1,:), coef_pol_x_init(1,:), coef_pol_dif_init(1,:), ...
    trans_init(1), list_amp_init(1,:), list_lambda, list_lambda_cal, ...
    rad_ROI, pix_IFS, lamp_spec, pattern_model, IFS_calib_wave, ...
    IFS_calib_spec, IFS_dark_wave, IFS_dark_spec, ...
    IFS_W, mu_wave_norm, option_opti_spec) ;

% Update of the simulation
IFS_sim_wave(list_i, list_j) = ...
    IFS_sim_wave(list_i, list_j) + Sim_wave ;
IFS_sim_spec(list_i, list_j) = ...
    IFS_sim_spec(list_i, list_j) + Sim_spec ;

ind_hex = pos2ind([coef_pol_y(1,1), coef_pol_x(1,1)], pix_map) ;
map_spec(ind_hex(1)+(-nb_rad*rad_ROI:nb_rad*rad_ROI), ...
    ind_hex(2)+(-nb_rad*rad_ROI:nb_rad*rad_ROI)) = 1 ;


% Estimation of the polynomial y-law
Pol_y = zeros(nb_lambda_cal, deg_pol_y+1) ;
for d = 0:deg_pol_y
    Pol_y(:,d+1) = list_lambda_cal.^d ;
end

% Estimation of the polynomial x-law
Pol_x = zeros(nb_lambda_cal, deg_pol_x+1) ;
for d = 0:deg_pol_x
    Pol_x(:,d+1) = list_lambda_cal.^d ;
end

%% Loop on the points to analyze
spec = 1  ;
nb_spec = 1 ;
ind_new = 1 ;

amp_spec = list_amp(spec,:) ;
coef_pol_y_spec = coef_pol_y(spec,:) ;
coef_pol_x_spec = coef_pol_x(spec,:) ;
coef_pol_dif_spec = coef_pol_dif(spec,:) ;

delta_per = 0.01 ;
[per_aux, lastprint] = display_percentage('init', ...
    'Fitting spectra') ;
while spec < nb_spec+1
    %% Display the percentage
    [per_aux, lastprint] = display_percentage('iter', ...
        {sum(map_spec(:)>0)/surf_hex*100, per_aux, delta_per, lastprint}) ;

    %% Finding the next hexagon to analyze
    % List of the local hexagon corners
    list_neigh = get_neigh_hex(map_spec, ...
        [coef_pol_y(spec,1), coef_pol_x(spec,1)], side_hex, theta_hex) ;

    % Neighboring corners
    list_c = (1:6)' ;
    list_c = list_c(list_neigh>0) ;
    list_c_new = setdiff((1:6)', list_c) ;


    %% Update of side_hex and theta_hex
    % Parameters of the local hexagon defined with all the correctly
    % positioned neighboring hexagons in all the wavelengths
    if ~isempty(list_c)
        list_neigh_c = list_neigh(list_c) ;
        % Averaging amplitudes in the haxagon
        amp_spec = mean(list_amp([list_neigh_c; spec],:), 1) ;
        coef_pol_y_spec = mean(coef_pol_y([list_neigh_c; spec],:), 1) ;
        coef_pol_x_spec = mean(coef_pol_x([list_neigh_c; spec],:), 1) ;
        coef_pol_dif_spec = mean(coef_pol_dif([list_neigh_c; spec],:), 1) ;
        
        % Refining the edge
        side_hex = (coef_pol_x(list_neigh_c,1)- coef_pol_x(spec,1)) + ...
            1i*(coef_pol_y(list_neigh_c,1)- coef_pol_y(spec,1)) ;
        
        side_hex = mean( abs(side_hex).* ...
            exp(1i*mod(angle(side_hex)*180/pi,60)*pi/180)) ;
        theta_hex_new = mod(180/pi*angle(side_hex)+30,60)-30 ;
            % Orientation of the hexagon (in [-30°, 30°])
        side_hex = abs(side_hex) ; % Side of the hexagon

        % Checking that the offset did not wrap of 60° because of value
        % close to -30° / 30°
        if abs(theta_hex_new-theta_hex)<30
            theta_hex = theta_hex_new ;
        end
    end
    
    % Current hexagon position
    coef_pol_y_spec(1) = coef_pol_y(spec,1) ;
    coef_pol_x_spec(1) = coef_pol_x(spec,1) ;
    
    
    %% Robust detection in all wavelength
    if ~isempty(list_c_new) && ind_new<length(list_c_new)+1
        %% Identification of the position of the new hexagon by fitting all
        % the wavelength at the same time
        
        % Creating the local hexagon
        c_new = list_c_new(ind_new) ;
        [x_c, y_c] = rot_2D(theta_hex+(c_new-1)*60, ...
            side_hex, 0) ;
        
        % Updating position
        coef_pol_y_spec(1) = y_c + coef_pol_y_spec(1) ;
        coef_pol_x_spec(1) = x_c + coef_pol_x_spec(1) ;

        inField = true ;
        flag_fit = true ;
        list_ind = pos2ind([coef_pol_y_spec(1), coef_pol_x_spec(1)], ...
            pix_IFS) ;
        
        % Checking on the background that the position is in the field
        % of view
        inField = inField & (list_ind(1)>=1) & ...
            (list_ind(1)<=pix_IFS.nb_y)  ; % y-axis
        inField = inField & (list_ind(2)>=1) & ...
            (list_ind(2)<=pix_IFS.nb_x)  ; % y-axis
        list_ind(1) = max(1, min(pix_IFS.nb_y, list_ind(1))) ;
        list_ind(2) = max(1, min(pix_IFS.nb_x, list_ind(2))) ;
        inField = inField & IFS_mask(list_ind(1), list_ind(2));
        
        if inField
            %% Position of the patterns before the fit
            list_y = Pol_y*coef_pol_y_spec' ;
            list_x = Pol_x*coef_pol_x_spec' ;
            
            %% Fitting spectrum
            [coef_pol_y_spec, coef_pol_x_spec, ...
                coef_pol_dif_spec, trans_spec, ...
                amp_spec, list_i, list_j, Sim_wave, Sim_spec] = ...
                fit_calib_spectrum( ...
                coef_pol_y_spec, coef_pol_x_spec, ...
                coef_pol_dif_spec, 1, ...
                amp_spec, list_lambda, list_lambda_cal, ...
                rad_ROI, pix_IFS, lamp_spec, ...
                pattern_model, IFS_calib_wave, IFS_calib_spec, ...
                IFS_dark_wave, IFS_dark_spec, ...
                IFS_W, mu_wave_norm, option_opti_spec) ;
            
            %% Checking if the results did not diverge beyond a rad_ROI
            % Position of the patterns
            list_y_out = Pol_y*coef_pol_y_spec' ;
            list_x_out = Pol_x*coef_pol_x_spec' ;
            delta_y = abs(list_y_out-list_y) ;
            delta_x = abs(list_x_out-list_x) ;
            
            flag_fit = tol_trans<trans_spec && max(delta_y)<rad_ROI && ...
                max(delta_x)<rad_ROI ;
        end

        %% Adding new hexagons in the list if needed
        if inField && flag_fit
            % Update of the number of  hexagons
            nb_spec = nb_spec + 1 ;

            % Update the polynomial laws
            coef_pol_y(nb_spec,:) = coef_pol_y_spec ;
            coef_pol_x(nb_spec,:) = coef_pol_x_spec ;
            coef_pol_dif(nb_spec,:) = coef_pol_dif_spec ;
            
            % Update the tranmission parameters
            trans(nb_spec) = trans_spec ;
            list_amp(nb_spec,:) = amp_spec ;

            % Update the index map
            ind_hex = pos2ind([coef_pol_y_spec(1), coef_pol_x_spec(1)], ...
                pix_map) ;
            
            % Update of the simulation
            IFS_sim_wave(list_i, list_j) = ...
                IFS_sim_wave(list_i, list_j) + Sim_wave ;
            IFS_sim_spec(list_i, list_j) = ...
                IFS_sim_spec(list_i, list_j) + Sim_spec ;

% % % % % % % % %                 % Checking there is no overwriting
            if map_spec(ind_hex(1), ind_hex(2))>0
                % Loop on the spectra
                list_y = zeros(nb_spec, nb_lambda_cal) ;
                list_x = zeros(nb_spec, nb_lambda_cal) ;
                for spec_aux = 1:nb_spec
                    % Position of the patterns
                    list_y(spec_aux,:) = Pol_y*coef_pol_y(spec_aux,:)' ;
                    list_x(spec_aux,:) = Pol_x*coef_pol_x(spec_aux,:)' ;
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
                title('Epic fail!') ;

                subplot(1,2,2) ;
                imshow(map_spec, []) ;
                colormap(gca,jet) ;
                title('Index map of the neighboring hexagons') ;
                
                figure(2) ;
                subplot(1,2,1) ;
                imshow(IFS_sim_wave, [val_min_wave, val_max_wave]) ;
                subplot(1,2,2) ;
                imshow(IFS_sim_spec, [val_min_spec, val_max_spec]) ;
                
                error('epic fail') ;
            end
% % % % % % % % %

            map_spec(...
                ind_hex(1)+(-nb_rad*rad_ROI:nb_rad*rad_ROI), ...
                ind_hex(2)+(-nb_rad*rad_ROI:nb_rad*rad_ROI)) = ...
                nb_spec ;
        else
            % Update of ind_new
            ind_new = ind_new+1 ;
        end
    else
        %% Update of hex and ind_new
        spec = spec+1 ;
        ind_new = 1 ;
    end
end
% Cleaning percentage
display_percentage('exit', lastprint) ;

% Keeping spectra
coef_pol_y = coef_pol_y(1:nb_spec,:) ;
coef_pol_x = coef_pol_x(1:nb_spec,:) ;
coef_pol_dif = coef_pol_dif(1:nb_spec,:) ;
trans = trans(1:nb_spec,:) ;
list_amp = list_amp(1:nb_spec,:) ;

%% Display current results
% Loop on the spectra
list_y = zeros(nb_spec, nb_lambda_cal) ;
list_x = zeros(nb_spec, nb_lambda_cal) ;
for spec_aux = 1:nb_spec
    % Position of the patterns
    list_y(spec_aux,:) = Pol_y*coef_pol_y(spec_aux,:)' ;
    list_x(spec_aux,:) = Pol_x*coef_pol_x(spec_aux,:)' ;
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

figure(2) ;
subplot(1,2,1) ;
imshow(IFS_sim_wave, [val_min_wave, val_max_wave]) ;
subplot(1,2,2) ;
imshow(IFS_sim_spec, [val_min_spec, val_max_spec]) ;


%% Saving current results
save([save_path, 'coef_pol_y.txt'], 'coef_pol_y', '-ascii', '-double') ;
save([save_path, 'coef_pol_x.txt'], 'coef_pol_x', '-ascii', '-double') ;
save([save_path, 'coef_pol_dif.txt'], 'coef_pol_dif', '-ascii', ...
    '-double') ;
save([save_path, 'trans.txt'], 'trans', '-ascii', '-double') ;
save([save_path, 'list_amp.txt'], 'list_amp', '-ascii', '-double') ;
save_fits(map_spec, 'map_spec', save_path) ;

%% IFS simulation
save_fits(IFS_sim_wave, 'IFS_sim_wave', save_path) ;
save_fits(IFS_sim_spec, 'IFS_sim_spec', save_path) ;
save_fits(IFS_calib_wave - IFS_sim_wave, 'IFS_res_wave', save_path) ;
save_fits(IFS_calib_spec - IFS_sim_spec, 'IFS_res_spec', save_path) ;