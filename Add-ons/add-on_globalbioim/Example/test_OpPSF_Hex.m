% This macro is a procedure to test the OpPSF_Hex operator
%
% This file supposes that the following toolboxes are in the Matlab path:
%   -> GlobalBioIm
%   -> add-on_GlobalBioIm
%
% Created: 02/22/2019 (mm/dd/yyyy)
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%


clear all ;
close all ;
clc ;

%% Test apply_ and downsampling
nb_p = 10 ;
rate_OS = 8 ;
x = (-nb_p:nb_p)' ;
y = (-nb_p:nb_p)' ;
y_c = 0.001 ;
x_c = 0*-5 ;
theta = 40 ;
sig = 3 ;
amp = 1 ;
offset = 0 ;
eps = 1e-32 ;

[OpDS, y_OS, x_OS] = get_downsampling_operator( ...
    rate_OS, y, x) ;
PSF_Hex = OpPSF_Hex(y_OS, x_OS, [], [], [], [], [], [], eps) ;
I = PSF_Hex*[y_c ; x_c ; theta ; sig ; amp ; offset] ;

figure(1)
subplot(1,2,1) ;
imshow(I, [offset, amp/sig^2/10]) ;
subplot(1,2,2) ;
imshow(OpDS*I, [offset, amp/sig^2/2]) ;
% imshow(OpDS*I, []) ;

disp(sum(I(:))) ;

%%
% pix.dx = 1 ;
% pix.dy = 1 ;
% pix.nb_x_pad = length(x) ;
% pix.nb_y_pad = length(y) ;
% flag_space = 'SD' ;
% const = 50   ;
% 
% [~, FT_hexagon] = get_conv_hexagon(pix, sqrt(2*pi)*nb_p/sig, theta, ...
%     1, flag_space, const) ;
% FT_hexagon = abs(FT_hexagon).^2 ;
% FT_hexagon = FT_hexagon./sum(FT_hexagon(:)) ;
% FT_hexagon = offset+amp*FT_hexagon ;
% 
% figure(1)
% subplot(1,2,1) ;
% imshow(I, [offset, amp/sig^2/100]) ;
% subplot(1,2,2) ;
% imshow(FT_hexagon, [offset, amp/sig^2/100]) ;
% 
% figure(2)
% imshow(I-FT_hexagon, [-amp/sig^2/100, amp/sig^2/100]) ;


%% Test optimization
nb_x = 128 ;
nb_y = 64 ;
x = get_Fourier_vector(nb_x, 1) ;
y = get_Fourier_vector(nb_y, 1) ;
y_c = 15.6 ;
x_c = 5.2 ;
theta = 40 ;
sig = 15 ;
amp = 1 ;
offset = 0.5 ;
eps = 1e-32 ;


% Differences for the initialization of the optimization
delta_0 = [-3, 5, 5, -2, -0.75, 3]' ; 
sig_SNR = 0.0001 ;
option_opti.method = 'VMLMB' ;
option_opti.verbose = true ;
option_opti.maxiter = 5000 ;

% Simulation
PSF_Hex = OpPSF_Hex(y, x, [], [], [], [], [], [], eps) ;

par_th = [y_c ; x_c ; theta ; sig ; amp ; offset] ;
I = PSF_Hex*par_th ;
I_s = I+sig_SNR*randn(size(I)) ;


checkMap(PSF_Hex) ;


% Test optimization
C = CostL2(PSF_Hex.sizeout,I_s)*PSF_Hex ;
const_min = [] ;
const_max = [] ;

[par_sol] = run_Opti(C, const_min, const_max, ...
    option_opti, par_th+delta_0) ;

disp('Theoritical values:') ;
disp(par_th) ;
disp('Estimated values:') ;
disp(par_sol) ;

figure(1)
subplot(2,2,1)
imshow(I, [offset, offset+amp/sig^2/10]) ;
title('Initial pattern') ;
subplot(2,2,2)
imshow(I_s, [offset, offset+amp/sig^2/10]) ;
title('Noisy pattern') ;
subplot(2,2,3)
imshow(PSF_Hex*par_sol, [offset, offset+amp/sig^2/10]) ;
title('Retrieved pattern') ;
subplot(2,2,4)
imshow(PSF_Hex*(par_sol+delta_0), ...
    [offset+delta_0(6), offset+delta_0(6)+amp/sig^2/10]) ;
title('Initial pattern') ;


%% Test optimization Gaussian
nb_x = 100 ;
nb_y = 32 ;
rate_OS = 4 ;
x = get_Fourier_vector(nb_x, 1)' ;
y = get_Fourier_vector(nb_y, 1)' ;
y_c = -3.6 ;
x_c = 2.2 ;
theta = 19 ;
sig = 4.5 ;
amp = 1 ;
offset = 0*0.5 ;
eps = 1e-32 ;

% Differences for the initialization of the optimization
delta_0 = 0*[-1, 0.5, 5, -2, -0.75, 3]' ; 
sig_SNR = 0.0001 ;
option_opti.method = 'VMLMB' ;
option_opti.verbose = true ;
option_opti.maxiter = 5000 ;


% Simulation
if rate_OS>1
    [OpDS, y_OS, x_OS] = get_downsampling_operator( ...
        rate_OS, y, x) ;
    PSF_Hex = OpDS*OpPSF_Hex(y_OS, x_OS, [], [], [], [], [], [], eps) ;
else
    PSF_Hex = OpPSF_Hex(y, x, [], [], [], [], [], [], eps) ;
end
% Gauss = OpDS*OpPSF_Hex(y_OS, x_OS, [], [], theta, [], [], [], eps) ;
% Gauss = OpDS*OpGauss(y_OS, x_OS, [], [], [], [], [], eps) ;
Gauss = OpGauss(y, x, [], [], [], [], [], true) ;



par_th = [y_c ; x_c ; theta ; sig ; amp ; offset] ;
I = PSF_Hex*par_th ;
I_s = I+sig_SNR*randn(size(I)) ;


% Test optimization
C = CostL2(PSF_Hex.sizeout,I_s)*Gauss ;
const_min = [] ;
const_max = [] ;

par_in = par_th+delta_0 ;
par_in = par_in([1,2,4,5,6]) ;
[par_sol] = run_Opti(C, const_min, const_max, ...
    option_opti, par_in) ;

disp('Theoritical values:') ;
disp(par_th) ;
disp('Estimated values:') ;
disp(par_sol) ;

figure(1)
subplot(2,2,1)
imshow(I, [offset, offset+amp/sig^2/5]) ;
title('Initial pattern') ;
subplot(2,2,2)
imshow(I_s, [offset, offset+amp/sig^2/5]) ;
title('Noisy pattern') ;
subplot(2,2,3)
imshow(Gauss*par_sol, [offset, offset+amp/sig^2/5]) ;
title('Retrieved pattern') ;
subplot(2,2,4)
imshow(I_s-Gauss*par_sol, ...
    [-amp/sig^2/10, amp/sig^2/5]) ;
title('Residues') ;

disp(par_th(4)/par_sol(3))