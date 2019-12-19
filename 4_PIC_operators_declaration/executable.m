%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This macro is a procedure to simulate the propagation model of SPHERE for
% a given experiment
%
% Created: 03/21/2018 (mm/dd/yyyy)
% Created: 08/08/2018 (mm/dd/yyyy) New functions
% Created: 11/16/2018 (mm/dd/yyyy) Adaptation to global laws
% Created: 11/27/2018 (mm/dd/yyyy) Adaptation to polynomial laws
% Created: 03/13/2019 (mm/dd/yyyy) Adaptation to local polynomial laws and
% autocalibration on the data
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
addpath(genpath('./functions')) ;
addpath('./procedures') ;

% Loading toolboxes
addpath(genpath('../Add-ons/')) ;


%% Set parameters
run parameters
save_path = [save_path, ...
    datestr(now,'yyyy_mm_dd_HH_MM_SS'),'/'] ;
save_path = make_dir(save_path) ;
save_path_coef_model = [save_path, 'coef_model/'] ;
if ~exist(save_path_coef_model, 'dir')
    save_path_coef_model = make_dir(save_path_coef_model) ;
end

% Saving the running codes
save_path_aux = make_dir([save_path, 'code/']) ;
copyfile('./executable.m', save_path_aux) ;
save_path_aux = make_dir([save_path, 'code/procedures/']) ;
copyfile('./procedures',save_path_aux) ;
save_path_aux = make_dir([save_path, 'code/functions/']) ;
copyfile('./functions',save_path_aux) ;


%% Initialization
run initialization


%% Estimating the shift
if flag_estimate_shift
    run estimate_shift
else
    delta_pos = [0;0] ;
end
coef_pol_x(:, 1) = coef_pol_x(:, 1)+delta_pos(2) ;
coef_pol_y(:, 1) = coef_pol_y(:, 1)+delta_pos(1) ;


%% Update the parameters
spectra.coef_pol_y = coef_pol_y ;
spectra.coef_pol_x = coef_pol_x ;


%% Saving model
% Coefficients of the polynomial fittings
save([save_path_coef_model, 'coef_pol_dif.txt'], 'coef_pol_dif', ...
    '-ascii', '-double') ;
save([save_path_coef_model, 'coef_pol_x.txt'], 'coef_pol_x', '-ascii', ...
    '-double') ;
save([save_path_coef_model, 'coef_pol_y.txt'], 'coef_pol_y', '-ascii', ...
    '-double') ;

% Wavelength
save([save_path, 'list_lambda.txt'], 'list_lambda', '-ascii', '-double') ;
save([save_path, 'lambda_0.txt'], 'lambda_0', '-ascii', '-double') ;

% Transmission
save([save_path, 'trans.txt'], 'trans', '-ascii', '-double') ;
save([save_path, 'lamp_spec.txt'], 'lamp_spec', '-ascii', '-double') ;


%% Model construction
disp('Model construction...') ;
tic_aux = tic ;

% Padding operator
fprintf('    Padding operator...') ;
tic
Pad_2D = LinOpPad([pix.nb_y, pix.nb_x], [pix.nb_y_pad, pix.nb_x_pad]) ;
save([save_path, 'Pad_2D.mat'], 'Pad_2D') ;
Pad = LinOpPad([pix.nb_y, pix.nb_x, nb_lambda], [pix.nb_y_pad, ...
    pix.nb_x_pad, nb_lambda]) ;
save([save_path, 'Pad.mat'], 'Pad') ;
disp(['   ', num2str(toc), ' s']) ;


% Convolution kernel
fprintf('    Convolution kernel...') ;
tic
Conv_hexa_2D = get_conv_hexagon(pix, lenslet_side, lenslet_theta, 1, ...
    flag_space, const, true) ;
save([save_path, 'Conv_hexa_2D.mat'], 'Conv_hexa_2D') ;
Conv_hexa = get_conv_hexagon(pix, lenslet_side, lenslet_theta, ...
    nb_lambda, flag_space, const, true) ;
save([save_path, 'Conv_hexa.mat'], 'Conv_hexa') ;
disp(['   ', num2str(toc), ' s']) ;


% Interpolator (positions in the sky assumed to be at the position of the
% absolute grid)
fprintf('    Interpolator...') ;
tic
Interp_2D = LinOpInterp({get_Fourier_vector(pix.nb_y, 1)', ...
    get_Fourier_vector(pix.nb_x, 1)'}, ...
    [list_y/pix.dy, list_x/pix.dx]) ;
save([save_path, 'Interp_2D.mat'], 'Interp_2D') ;
Interp = LinOpInterp({get_Fourier_vector(pix.nb_y, 1)', ...
    get_Fourier_vector(pix.nb_x, 1)'}, ...
    [list_y/pix.dy, list_x/pix.dx], [], ...
    [pix.nb_y, pix.nb_x, nb_lambda], [1,2]) ;
save([save_path, 'Interp.mat'], 'Interp') ;
disp(['   ', num2str(toc), ' s']) ;

% Sparse matrix projector
fprintf('    Sparse matrix projector...') ;
tic
Spec_proj = LinOpSpecProjLocLaw(sensor, spectra, pattern_model) ;
disp(['   ', num2str(toc), ' s']) ;
fprintf('    Sparse matrix projector saving...') ;
tic
try
    save([save_path, 'Spec_proj.mat'], 'Spec_proj') ;
catch
    disp('   Saving in -v7.3 format') ;
    save([save_path, 'Spec_proj.mat'], 'Spec_proj','-v7.3') ;
end
disp(['   ', num2str(toc), ' s']) ;


fprintf('Initialization done!') ;
disp(['   ', num2str(toc(tic_aux)), ' s']) ;

%% Plot of the convolution kernel
nb_h_plot = 37 ;
fig = figure(1) ;
test = zeros([pix.nb_y_pad, pix.nb_x_pad]) ;
test(1,1) = 1 ;
test = fftshift(test) ;
pix_test = pix ;
pix_test.nb_x = pix_test.nb_x_pad ;
pix_test.nb_y = pix_test.nb_y_pad ;
imshow(Conv_hexa_2D*test, []) ;
plot_hexagon(fig, ...
    [list_y(1:nb_h_plot)-list_y(1:nb_h_plot), ...
    list_x(1:nb_h_plot)-list_x(1:nb_h_plot)], ...
    lenslet_theta, lenslet_side, pix_test, 'r', LineWidth) ;

%% Simulation
hypercube = repmat(reshape(lamp_spec, [1,1,nb_lambda]), ...
    [pix.nb_y, pix.nb_x, 1]) ;

disp('Simulation...') ;
tic_aux = tic ;
% Convolution
fprintf('   Convolution...') ;
tic
hypercube_conv = Pad'*Conv_hexa*Pad*hypercube ;
disp(['   ', num2str(toc), ' s']) ;

% Interpolation
fprintf('   Interpolation...') ;
tic
hypercube_interp = Interp*hypercube_conv ;
disp(['   ', num2str(toc), ' s']) ;

% Projection
fprintf('   Projection...') ;
tic
data_sim = Spec_proj*hypercube_interp ;
disp(['   ', num2str(toc), ' s']) ;

fig = figure(2) ;
% imshow(abs(data_sim).^0.25, []) ;
imshow(data_sim, []) ;
plot_dot(fig, ...
    [list_y, ...
    list_x], ...
    sensor, [1,0.5,0.5], MarkerSize) ;

plot_hexagon(fig, ...
    [list_y(1:nb_h_plot), ...
    list_x(1:nb_h_plot)], ...
    lenslet_theta, lenslet_side, sensor, 'w', LineWidth) ;

fprintf('Simulation done!') ;
disp(['   ', num2str(toc(tic_aux)), ' s']) ;

%% Saving
save_fits(data_sim, 'data_sim', save_path) ;