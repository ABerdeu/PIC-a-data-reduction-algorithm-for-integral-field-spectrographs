%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This macro is a procedure to analyze the flat calibration
%
% Created: 09/07/2018 (mm/dd/yyyy)
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all ;
clc ;

% Absolute path
abs_path = pwd ;
abs_path = [abs_path, '\'] ;

% Load functions and packages
restoredefaultpath

% Loading toolboxes
addpath(genpath('../Add-ons/')) ;


%% Set parameters
% Calibration files
path_cal = './results/path_to_the_reduction/' ;
nb_fit = 10 ;
dyn_res = 1000 ;
dyn_res_white = 20 ;
nb_bin = 256 ;

%% Opening data
path_ref = [path_cal, 'refinement_pos_dif_par/'] ;
suf_fit = ['_it_', num2str(nb_fit)] ;

list_y = load([path_ref, ['coef/coef_pol_y', suf_fit, '.txt']], '-ascii') ;
list_x = load([path_ref, ['coef/coef_pol_x', suf_fit, '.txt']], '-ascii') ;

list_y = list_y(:,1) ;
list_x = list_x(:,1) ;

global lamp_spec ;       
global trans ;      

% Parameters
lamp_spec = load([path_cal, 'lamp_spec.txt'], '-ascii') ;
trans = load([path_cal, 'trans.txt'], '-ascii') ;
nb_lenslet = length(trans) ;
nb_lambda = length(lamp_spec) ;

list_lamp_spec = repmat(lamp_spec, [1, nb_fit+1]) ;
list_trans = repmat(trans, [1, nb_fit+1]) ;
for it = 1:nb_fit
    % Loading f^th fit
    lamp_spec = load([path_ref, 'coef/lamp_spec_it_', num2str(it), ...
        '.txt'], '-ascii') ;
    trans = load([path_ref, 'coef/trans_it_', num2str(it), '.txt'], ...
        '-ascii') ;
    list_lamp_spec(:,it+1) = lamp_spec ;
    list_trans(:,it+1) = trans ;
end


%% Analyze spectral response
% Plot parameters
LineWidth = 3 ;
fig_pos = [75, 100] ;
fig_size = 850*[1.5, 1] ;
list_color = jet(nb_fit+1) ;
FontSize = 30 ;
FontSize_axis = 20 ;

% Display
fig_1 = figure(1) ;
set(fig_1,'rend','painters','pos', [fig_pos, fig_size]) ;
clf(fig_1) ;
hold on
for l = 1:nb_fit+1
    plot(1:nb_lambda, list_lamp_spec(:,l), '-', 'Color', ...
        list_color(l,:), 'LineWidth', LineWidth) ;
end
hold off
axis([1, nb_lambda, 0, 4*10^4]) ;
set(gca,'FontSize', FontSize_axis);
xlabel('$\lambda_i$', 'FontSize', FontSize) ;
ylabel('adu', 'FontSize', FontSize) ;
title('Spectral response along iterations', 'FontSize', FontSize) ;
set(findall(gcf,'-property','Interpreter'),'Interpreter','Latex') ;

% Saving
saveas(fig_1, [path_ref, 'Spectral_fit'], 'png') ;

%% Analyze transmission
% Parameters
MarkerSize = 12 ; % Marker size
vect_view = [-1 -3 1] ; % Viewing vector
nb_color = 1024 ;
z_axis = 1 + 0.15*[-1, 1] ;
color_axis = 1 + [-0.075, 0.075] ;
fig_pos = [75, 100] ;
fig_size = 900*[2, 1] ;
FontSize_axis = 15 ;
FontSize = 30 ;
ind_scale_bar = 150:1200 ; % Index to insert the scale bar

option_scalebar.nb_pix_bar = 49 ;
option_scalebar.nb_pix_delta = 0 ;
option_scalebar.colormap = jet(1024) ;
option_scalebar.nb_stick = 10 ;
option_scalebar.str_format = '%0.3e' ;
option_scalebar.FontSize = 40 ;
option_scalebar.BGcolor = [1,1,1] ;
option_scalebar.nb_edge = 0 ;

% Plotting options
min_axis = min(color_axis) ;
max_axis = max(color_axis) ;
list_color = jet(nb_color) ;
ind_color = list_trans(:,nb_fit+1) ;
ind_color = round((ind_color-min_axis) / ...
    (max_axis-min_axis)*nb_color) ;
ind_color = min(max(1, ind_color), nb_color) ;
list_color = list_color(ind_color, :) ;

% Display
fig_2 = figure(2) ;
clf(fig_2) ;
set(fig_2,'rend','painters','pos', [fig_pos, [fig_size(2), fig_size(2)]]) ;
scatter(list_x, list_y, MarkerSize, list_color, 'filled') ;
axis equal ;
axis([-1024, 1024, -1024, 1024]) ;
set(gca,'FontSize',FontSize_axis);
title('Transmission for each lenslet', 'FontSize', FontSize) ;
set(findall(gca,'-property','Interpreter'),'Interpreter','Latex') ;
set(findall(gcf,'-property','Interpreter'),'Interpreter','Latex') ;
view([0,0,1]) ;
saveas(fig_2, [path_ref, 'Transmission_fit'], 'png') ;

% Inserting scale bar
pic_save = imread([path_ref, 'Transmission_fit.png']) ;
pic_scale = insert_colorbar_pic( ...
    linspace(max_axis, min_axis, length(ind_scale_bar))', [], ...
    option_scalebar) ;
pic_save = cat(2, pic_save, 255*ones(size(pic_save,1), ...
    size(pic_scale,2), 3)) ;
delta_ind = size(pic_scale,1)-length(ind_scale_bar) ;
delta_ind = (-delta_ind+round(delta_ind/2) + ind_scale_bar(1)): ...
    (ind_scale_bar(end)+delta_ind-round(delta_ind/2)+1) ;
pic_save(delta_ind, (end-size(pic_scale,2)+1):end,:) = ...
    255*pic_scale ;

% Saving
imwrite(pic_save, [path_ref, 'Transmission_fit'], 'png') ;

%% Analayzing the residues statistics
disp('Analyzing residues...') ;
run([path_cal, 'code/procedures/parameters.m']) ;
list_color = jet(nb_fit+1) ;

% Loading variables
IFS_calib_spec = fitsread([path_cal, 'IFS_calib_spec_corrected.fits']) ;
IFS_mask = fitsread([data.sensor_calib, 'IFS_mask.fits']) ;
IFS_BP = fitsread([path_cal, 'IFS_BP_corrected.fits']) ;
IFS_sensor_flat = fitsread([data.sensor_calib, 'IFS_sensor_flat.fits']) ;
IFS_dark_current = fitsread([data.sensor_calib, 'IFS_dark_BG.fits']) ;
IFS_dark_spec = get_fits_exposure(data.calib_spec)*IFS_dark_current ./ ...
    IFS_sensor_flat ;


% Selecting the pixels
pix_sel = IFS_mask & IFS_BP ;

% Bins
list_edges = dyn_res*linspace(-1, 1, nb_bin) ;
list_edges_white = dyn_res_white*linspace(-1, 1, nb_bin) ;
list_bin = (list_edges(2:nb_bin) + list_edges(1:(nb_bin-1)))/2 ;
list_bin_white = (list_edges_white(2:nb_bin) + ...
    list_edges_white(1:(nb_bin-1)))/2 ; 

% Noise law
if flag_MAD
    eta = load([data.sensor_calib, 'eta_MAD.txt'], '-ascii') ;
    var_0 = load([data.sensor_calib, 'sig_MAD.txt'], '-ascii') ;
else
    eta = load([data.sensor_calib, 'eta.txt'], '-ascii') ;
    var_0 = load([data.sensor_calib, 'sig.txt'], '-ascii') ;
end

% Building histograms
N = zeros(nb_bin-1, nb_fit+1) ;
sig = zeros(nb_fit+1, 1) ;
N_white = zeros(nb_bin-1, nb_fit+1) ;
sig_white = zeros(nb_fit+1, 1) ;
for it = 0:nb_fit
    disp(['   Iteration: ', num2str(it), '/', num2str(nb_fit)]) ;
    
    % Loading simulated signal
    IFS_sim_spec = fitsread([path_ref, 'spec/IFS_sim_spec_it_', ...
        num2str(it), '.fits']) ;
    
    % Residues
    IFS_res_spec = IFS_calib_spec-IFS_sim_spec ;
    
    % Whitening the residues with the Poisson noise model
    IFS_res_spec_white = IFS_res_spec ./ ...
        (eta*(IFS_sim_spec+IFS_dark_spec)+var_0).^0.5 ;
    
    % Counting
    sig_aux = IFS_res_spec(pix_sel) ;
    sig_aux_white = IFS_res_spec_white(pix_sel) ;
    N(:,it+1) = histcounts(sig_aux, list_edges) ;
    N_white(:,it+1) = histcounts(sig_aux_white, list_edges_white) ;
    
    % Standard deviation
    sig_aux = sig_aux((sig_aux>-dyn_res)&(sig_aux<+dyn_res)) ;
    sig(it+1) = std(sig_aux(:)) ;
    sig_aux_white = ...
        sig_aux_white( ...
        (sig_aux_white>-dyn_res_white)&(sig_aux_white<+dyn_res_white)) ;
    sig_white(it+1) = std(sig_aux_white(:)) ;
end

% Display
fig_3 = figure(3) ;
set(fig_3,'rend','painters','pos', [fig_pos, fig_size]) ;
clf(fig_3) ;
plot(list_bin, ...
    max(N(:))*exp(-list_bin.^2/(2*sig(nb_fit+1)^2)) , ...
    '-', 'Color', [0,0,0], ...
    'LineWidth', LineWidth) ;
hold on
for it = 1:nb_fit+1
    plot(list_bin, N(:,it), '-', 'Color', list_color(it,:), ...
        'LineWidth', LineWidth) ;
end
hold off
axis([-dyn_res, +dyn_res, 0, max(N(:))]) ;
set(gca,'FontSize', FontSize_axis);
xlabel('$adu$', 'FontSize', FontSize) ;
ylabel('counts', 'FontSize', FontSize) ;
title('Residues histogram along iterations', 'FontSize', FontSize) ;
set(findall(gcf,'-property','Interpreter'),'Interpreter','Latex') ;

% Saving
saveas(fig_3, [path_ref, 'Residues_histogram'], 'png') ;

% Display
fig_4 = figure(4) ;
set(fig_4,'rend','painters','pos', [fig_pos, fig_size]) ;
clf(fig_4) ;
plot(list_bin_white, ...
    max(N_white(:))*exp(-list_bin_white.^2/(2*3^2)) , ...
    '-', 'Color', [0,0,0], ...
    'LineWidth', LineWidth) ;
hold on
for it = 1:nb_fit+1
    plot(list_bin_white, N_white(:,it), '-', 'Color', list_color(it,:), ...
        'LineWidth', LineWidth) ;
end
hold off
axis([-dyn_res_white, +dyn_res_white, 0, max(N_white(:))]) ;
set(gca,'FontSize', FontSize_axis);
xlabel('$adu$', 'FontSize', FontSize) ;
ylabel('counts', 'FontSize', FontSize) ;
title('Whitened residues histogram along iterations', 'FontSize', ...
    FontSize) ;
set(findall(gcf,'-property','Interpreter'),'Interpreter','Latex') ;

% Saving
saveas(fig_4, [path_ref, 'Residues_histogram_white'], 'png') ;