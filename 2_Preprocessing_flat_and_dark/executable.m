%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This macro is a procedure to analyze the flat and dark acquisitions
%
% Created: 01/26/2019 (mm/dd/yyyy)
% Modified: 02/25/2019 (mm/dd/yyyy) Analyzing all dark current
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
addpath('./procedures') ;
addpath('./functions') ;

% Loading toolboxes
addpath(genpath('../Add-ons/')) ;

%% Set parameters
run parameters
save_path = [save_path, ...
    datestr(now,'yyyy_mm_dd_HH_MM_SS'),'/'] ;
make_dir(save_path) ;
if save_lin_fits
    make_dir([save_path, 'dark/']) ;
    make_dir([save_path, 'dark_BG/']) ;
    make_dir([save_path, 'dark_flat/']) ;
    make_dir([save_path, 'sky_BG/']) ;
end
make_dir([save_path, 'figures/']) ;

% Saving the running codes
save_path_aux = make_dir([save_path, 'code/']) ;
copyfile('./executable.m', save_path_aux) ;
save_path_aux = make_dir([save_path, 'code/procedures/']) ;
copyfile('./procedures',save_path_aux) ;
save_path_aux = make_dir([save_path, 'code/functions/']) ;
copyfile('./functions',save_path_aux) ;

%% Analyzing spectral flat
run Analyze_spectral_flat

%% Analyzing darks
run Analyze_dark

%% Analyzing darks
run Analyze_dark_BG

%% Analyzing flats
run Analyze_flat

%% Global analysis to find the bad pixels
BP_mask = BP_mask_dark.*BP_mask_dark_BG.*BP_mask_flat ;

%% Analyzing sky backgrounds
run Analyze_sky_BG

%% Removing edges
BP_mask(1:nb_BP_edges, :) = 0 ;
BP_mask(:, 1:nb_BP_edges) = 0 ;
BP_mask(pix_IFS.nb_y+1+(-nb_BP_edges:-1), :) = 0 ;
BP_mask(:, pix_IFS.nb_x+1+(-nb_BP_edges:-1)) = 0 ;

save_fits(uint8(BP_mask), 'IFS_BP_mask', save_path) ;

%% Global analysis to find the noise law on the darks
disp('Analyzing noise law on the darks...') ;

% Dark
sig_dark_med = median(sig_dark(:)) ;
eta_dark_med = median(eta_dark(:)) ;
[list_med_dark, list_var_med_dark, list_var_MAD_dark] = ...
    get_variance_estimate(list_avg_dark, list_var_dark) ;

[fig_1, fig_2] = ...
    plot_linear_law(list_exp_time_dark, ...
    list_med_dark, eta_dark_med*list_exp_time_dark+sig_dark_med, ...
    list_var_med_dark, ...
    list_var_MAD_dark, 1, fig_pos, fig_size, LineWidth, MarkerSize, ...
    FontSize_axis, FontSize, [1,0,0]) ;
saveas(fig_1, [save_path, 'figures/fit_dark_adu.png']) ;
saveas(fig_2, [save_path, 'figures/fit_dark_exp_temp.png']) ;


%% Global analysis to find the noise law on the background darks
disp('Analyzing noise law on the background darks...') ;

% Dark
sig_dark_BG_med = median(sig_dark_BG(:)) ;
eta_dark_BG_med = median(eta_dark_BG(:)) ;
[list_med_dark_BG, list_var_med_dark_BG, list_var_MAD_dark_BG] = ...
    get_variance_estimate(list_avg_dark_BG, list_var_dark_BG) ;

[fig_3, fig_4] = ...
    plot_linear_law(list_exp_time_dark_BG, ...
    list_med_dark_BG, ...
    eta_dark_BG_med*list_exp_time_dark_BG+sig_dark_BG_med, ...
    list_var_med_dark_BG, ...
    list_var_MAD_dark_BG, 3, fig_pos, fig_size, LineWidth, MarkerSize, ...
    FontSize_axis, FontSize, [0,1,0]) ;

saveas(fig_3, [save_path, 'figures/fit_dark_BG_adu.png']) ;
saveas(fig_4, [save_path, 'figures/fit_dark_BG_exp_temp.png']) ;

%% Global analysis to find the noise law on the flats
disp('Analyzing noise law on the flats...') ;
% Flat
sig_flat_med = median(sig_flat(:)) ;
eta_flat_med = median(eta_flat(:)) ;
[list_med_flat, list_var_med_flat, list_var_MAD_flat] = ...
    get_variance_estimate(list_avg_flat, list_var_flat) ;

[fig_5, fig_6] = ...
    plot_linear_law(list_exp_time_flat, ...
    list_med_flat, ...
    eta_flat_med*list_med_flat+sig_flat_med, ...
    list_var_med_flat, ...
    list_var_MAD_flat, 5, fig_pos, fig_size, LineWidth, MarkerSize, ...
    FontSize_axis, FontSize, [0,0,1]) ;

saveas(fig_5, [save_path, 'figures/fit_flat_adu.png']) ;
saveas(fig_6, [save_path, 'figures/fit_flat_exp_temp.png']) ;


%% Fitting a global law
% Concatening data
list_exp_time = [ ...
    list_exp_time_dark; ...
    list_exp_time_dark_BG; ...
    list_exp_time_flat] ;

list_med = [ ...
    list_med_dark; ...
    list_med_dark_BG; ...
    list_med_flat] ;

list_var_med = [ ...
    list_var_med_dark; ...
    list_var_med_dark_BG; ...
    list_var_med_flat] ;

list_var_MAD = [ ...
    list_var_MAD_dark; ...
    list_var_MAD_dark_BG; ...
    list_var_MAD_flat] ;

% Reshaping
nb_acqui = length(list_exp_time) ;
list_exp_time = reshape(list_exp_time, [1, 1, nb_acqui]) ;
list_med = reshape(list_med, [1, 1, nb_acqui]) ;
list_var_med = reshape(list_var_med, [1, 1, nb_acqui]) ;
list_var_MAD = reshape(list_var_MAD, [1, 1, nb_acqui]) ;

% Fitting law
[a_var_MAD, b_var_MAD] = ...
    fit_linear_law(list_med, list_var_med, 1./list_var_MAD.^2, ...
    true) ;
[a_var, b_var] = ...
    fit_linear_law(list_med, list_var_med, [], ...
    true) ;

fig_7 = figure(7) ;
clf(fig_7) ;
set(fig_7,'rend','painters','pos', [fig_pos, fig_size]) ;
hold on
plot(list_med(:), a_var_MAD*list_med(:)+b_var_MAD, ...
    'LineWidth', LineWidth, 'MarkerSize', MarkerSize, ...
    'Color', [0, 0, 0]) ;
plot(list_med(:), a_var*list_med(:)+b_var, ...
    '--', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize, ...
    'Color', [0, 0, 0]) ;
errorbar(list_med_dark, list_var_med_dark, list_var_MAD_dark, ...
    'o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize, ...
    'Color', [1,0,0]) ;
errorbar(list_med_dark_BG, list_var_med_dark_BG, list_var_MAD_dark_BG, ...
    'o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize, ...
    'Color', [0,1,0]) ;
errorbar(list_med_flat, list_var_med_flat, list_var_MAD_flat, ...
    'o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize, ...
    'Color', [0,0,1]) ;
hold off
legend({'Weighted estimation', 'Unweighted estimation', ...
    'Dark (close)', 'Dark (open)', 'Flat'}, 'Location', 'northwest') ;
axis([0, 10*floor(max(list_med(:))/10+1), ...
    0, max(list_var_med_flat+list_var_MAD_flat)]) ;
set(gca,'FontSize',FontSize_axis);
set(findall(gca,'-property','Interpreter'),'Interpreter','Latex') ;
xlabel('Median flux (adu)', 'FontSize', FontSize) ;
ylabel('Variance (adu)', 'FontSize', FontSize) ;
set(findall(gcf,'-property','Interpreter'),'Interpreter','Latex') ;
saveas(fig_7, [save_path, 'figures/fit_combined_adu.png']) ;


% Saving fitted law
save([save_path, 'eta_MAD.txt'], 'a_var_MAD', '-ascii', '-double') ;
save([save_path, 'eta.txt'], 'a_var', '-ascii', '-double') ;
save([save_path, 'sig_MAD.txt'], 'b_var_MAD', '-ascii', '-double') ;
save([save_path, 'sig.txt'], 'b_var', '-ascii', '-double') ;

% Saving values to generate the figure
list_ext = {'_dark', '_dark_BG', '_flat'} ;
list_var = {'list_exp_time', 'list_med', 'list_var_MAD', 'list_var_med'} ;
for v = 1:length(list_var)
    for e = 1:length(list_ext)
        save([save_path, 'figures/', list_var{v}, list_ext{e}, '.txt'], ...
            [list_var{v}, list_ext{e}], '-ascii', '-double') ;
    end
end


%% Saving temporal darks and flats
save_fits(a_flat, 'IFS_sensor_flat', save_path) ;
save_fits(a_dark_BG, 'IFS_dark_BG', save_path) ;
save_fits(a_sky_BG, 'IFS_sky_BG', save_path) ;