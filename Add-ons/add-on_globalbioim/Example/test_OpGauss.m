% This macro is a procedure to test the OpGauss operator
%
% This file supposes that the following toolboxes are in the Matlab path:
%   -> GlobalBioIm
%   -> add-on_GlobalBioIm
%
% Created: 07/26/2018 (mm/dd/yyyy)
% Modified: 09/17/2018 (mm/dd/yyyy) Git desposit
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%


clear all ;
close all ;
clc ;

%% Parameters
x = (-5:5)' ;
y = (-10:10)' ;
y_c = 1 ;
x_c = 2 ;
theta = 30 ;
sig_par = 1.75 ;
sig_perp = 1.15 ;
amp = 4 ;
offset = 1 ;



%% On position, amplitude and offset
delta_0 = [0.75, 1, 1, -2]' ; % differences for the initialization of the
% optimization
sig_SNR = 0.1 ;
opti_method = 'Cauchy' ; % 'Cauchy', 'Huber', ...

% Simulation
Gauss = OpGauss(y, x, [], [], theta, sig_par, sig_perp, [], [], false) ;

I = Gauss*[y_c; x_c; amp; offset] ;
I_s = I+sig_SNR*randn(size(I)) ;
I_s=imnoise(I_s/(amp+offset), 'salt & pepper', 0.05) ;
I_s=I_s*(amp+offset) ;


checkMap(Gauss) ;


% Test optimization
s = 10 ;
% s = max(1.4826*median(abs(I_s(:)-median(I_s(:)))),0.002) ;
options.method = opti_method ;
options.flag_s = s ;
C = CostRobustPenalization(Gauss, I_s, options) ;
const_min = [-Inf; -Inf; 0.1; -Inf] ;
const_max = [+Inf; +Inf; +Inf; +Inf] ;
VMLMB=OptiVMLMB(C,const_min,const_max);  
VMLMB.ItUpOut=1; 
VMLMB.verbose=true ;
VMLMB.m=3;              % number of memorized step in hessian approximation
par_th = [y_c; x_c; amp; offset] ;
VMLMB.run(par_th+delta_0);

par_sol = VMLMB.xopt ;

disp('Theoritical values:') ;
disp(par_th) ;
disp('Estimated values:') ;
disp(par_sol) ;

figure(1)
subplot(2,2,1)
imshow(I, []) ;
title('Initial pattern') ;
subplot(2,2,2)
imshow(I_s, []) ;
title('Noisy pattern') ;
subplot(2,2,3)
imshow(Gauss*par_sol, []) ;
title('Retrieved pattern') ;
subplot(2,2,4)
imshow(C.computeW_(par_sol), []) ;
title('Equivalent weights') ;

%% On theta and sigma and amplitude and offset
delta_0 = [-2, 0.5, -0.25, 3, -2]' ; % differences for the initialization
% of the optimization
sig_SNR = 0.05 ;
opti_method = 'Cauchy' ; % 'Cauchy', 'Huber', ...

% Simulation
Gauss = OpGauss(y, x, y_c, x_c, [], [], [], [], []) ;

I = Gauss*[theta; sig_par; sig_perp; amp; offset] ;
I_s = I+sig_SNR*randn(size(I)) ;
I_s=imnoise(I_s/(amp+offset), 'salt & pepper', 0.05) ;
I_s=I_s*(amp+offset) ;

% Test optimization
s = 1.4826*median(abs(I_s(:)-median(I_s(:)))) ;
% s = 5 ;
options.method = opti_method ;
options.flag_s = s ;
C = CostRobustPenalization(Gauss, I_s, options) ;
const_min = [-90; 0.1; 0.1; 0; -Inf] ;
const_max = [+90; +Inf; +Inf; +Inf; +Inf] ;
VMLMB=OptiVMLMB(C,const_min,const_max);  
VMLMB.ItUpOut=1; 
VMLMB.verbose=false ;
VMLMB.maxiter=100000;      % max number of iterations
VMLMB.m=3;              % number of memorized step in hessian approximation
par_th = [theta; sig_par; sig_perp; amp; offset] ;
VMLMB.run(par_th+delta_0);

par_sol = VMLMB.xopt ;

disp('Theoritical values:') ;
disp(par_th) ;
disp('Estimated values:') ;
disp(par_sol) ;

figure(2)
subplot(2,2,1)
imshow(I, []) ;
title('Initial pattern') ;
subplot(2,2,2)
imshow(I_s, []) ;
title('Noisy pattern') ;
subplot(2,2,3)
imshow(Gauss*par_sol, []) ;
title('Retrieved pattern') ;
subplot(2,2,4)
imshow(C.computeW_(par_sol), []) ;
title('Equivalent weights') ;


%% On position, amplitude and elongation
delta_0 = [-0.25, 0.5, -0.25, 0.75, 3]' ; % differences for the initialization
% of the optimization
sig_SNR = 0.05 ;
opti_method = 'Cauchy' ; % 'Cauchy', 'Huber', ...

% Simulation
Gauss = OpGauss(y, x, [], [], theta, [], [], [], offset, false) ;

I = Gauss*[y_c; x_c; sig_par; sig_perp; amp] ;
I_s = I+sig_SNR*randn(size(I)) ;
I_s=imnoise(I_s/(amp+offset), 'salt & pepper', 0.1) ;
I_s=I_s*(amp+offset) ;


checkMap(Gauss) ;

% Test optimization
% s = max(1.4826*median(abs(I_s(:)-median(I_s(:)))),0.002) ;
s = 1 ;
options.method = opti_method ;
options.flag_s = s ;
C = CostRobustPenalization(Gauss, I_s, options) ;
const_min = [-Inf; -Inf; 0.1; 0.1; 0] ;
const_max = [+Inf; +Inf; +Inf; +Inf; +Inf] ;
VMLMB=OptiVMLMB(C,const_min,const_max);  
VMLMB.ItUpOut=1; 
VMLMB.verbose=false ;
VMLMB.maxiter=2000;      % max number of iterations
VMLMB.m=3;              % number of memorized step in hessian approximation
par_th = [y_c; x_c; sig_par; sig_perp; amp] ;
VMLMB.run(par_th+delta_0);

par_sol = VMLMB.xopt ;

disp('Theoritical values:') ;
disp(par_th) ;
disp('Estimated values:') ;
disp(par_sol) ;


figure(3)
subplot(2,2,1)
imshow(I, []) ;
title('Initial pattern') ;
subplot(2,2,2)
imshow(I_s, []) ;
title('Noisy pattern') ;
subplot(2,2,3)
imshow(Gauss*par_sol, []) ;
title('Retrieved pattern') ;
subplot(2,2,4)
imshow(C.computeW_(par_sol), []) ;
title('Equivalent weights') ;


%% On position, amplitude and elongation for axisymmetric
delta_0 = [-3, 0.5, 0.75, 3]' ; % differences for the initialization
% of the optimization
sig_SNR = 0.01 ;
opti_method = 'Cauchy' ; % 'Cauchy', 'Huber', 'Beaton-Tukey', ...

par_th = [y_c; x_c; sig_par; amp] ;
par_in = par_th+delta_0 ;

Gauss = OpGauss(y, x, [], [], [], [], offset) ;

I = Gauss*[y_c; x_c; sig_par; amp] ;
I_s = I+sig_SNR*randn(size(I)) ;
I_s=imnoise(I_s/(amp+offset), 'salt & pepper', 0.1) ;
I_s=I_s*(amp+offset) ;


checkMap(Gauss) ;

% Test optimization
% s = max(1.4826*median(abs(I_s(:)-median(I_s(:)))),0.002) ;
s = 1 ;
options.method = opti_method ;
options.flag_s = s ;
C = CostRobustPenalization(Gauss, I_s, options) ;
const_min = [] ;
const_max = [] ;
% const_min = [-Inf; -Inf; 0.1; 0] ;
% const_max = [+Inf; +Inf; +Inf; +Inf] ;
VMLMB=OptiVMLMB(C,const_min,const_max);  
VMLMB.ItUpOut=1; 
VMLMB.verbose=true ;
VMLMB.maxiter=10000;      % max number of iterations
VMLMB.m=3;              % number of memorized step in hessian approximation
VMLMB.run(par_in);

par_sol = VMLMB.xopt ;

disp('Theoritical values:') ;
disp(par_th) ;
disp('Estimated values:') ;
disp(par_sol) ;

figure(5)
subplot(2,2,1)
imshow(I, []) ;
title('Initial pattern') ;
subplot(2,2,2)
imshow(I_s, []) ;
title('Noisy pattern') ;
subplot(2,2,3)
imshow(Gauss*par_sol, []) ;
title('Retrieved pattern') ;
subplot(2,2,4)
imshow(C.computeW_(par_sol), []) ;
title('Equivalent weights') ;

%% Test Poisson noise (position, amplitude, elongation [axisymmetric])
delta_0 = 0.5*[-0.25, 0.5, 0.75, 3]' ; % differences for the initialization
% of the optimization
sig_SNR = 0.050 ;
eta_Poisson = 0.01 ;
opti_method = 'Cauchy' ; % 'Cauchy', 'Huber', 'Beaton-Tukey', ...

Gauss = OpGauss(y, x, [], [], [], [], offset) ;

I = Gauss*[y_c; x_c; sig_par; amp] ;
I_s = I+(sig_SNR+(eta_Poisson*(I-offset)).^0.5).*randn(size(I)) ;
I_s = I_s + (amp*imnoise(zeros(size(I_s)), 'salt & pepper', 0.05)) ;
checkMap(Gauss) ;

% Test optimization
options = [] ;
options.method = opti_method ;
% options.flag_s = 'none' ;
options.var_0 = sig_SNR^2-eta_Poisson*offset ;
options.eta = eta_Poisson ;
options.noise_model = 'Poisson' ;
C = CostRobustPenalization(Gauss, I_s, options) ;
const_min = [-Inf; -Inf; 0.1; 0] ;
const_max = [+Inf; +Inf; +Inf; +Inf] ;
VMLMB=OptiVMLMB(C,const_min,const_max);  
VMLMB.ItUpOut=1; 
VMLMB.verbose=true ;
VMLMB.maxiter=10000;      % max number of iterations
VMLMB.m=3;              % number of memorized step in hessian approximation
par_th = [y_c; x_c; sig_par; amp] ;
VMLMB.run(par_th+delta_0);

par_sol = VMLMB.xopt ;

disp('Theoritical values:') ;
disp(par_th) ;
disp('Estimated values:') ;
disp(par_sol) ;

figure(5)
subplot(2,2,1)
imshow(I, [min(I(:)), max(I(:))]+3*sig_SNR*[-1,1]) ;
title('Initial pattern') ;
subplot(2,2,2)
imshow(I_s, [min(I(:)), max(I(:))]+3*sig_SNR*[-1,1]) ;
title('Noisy pattern') ;
subplot(2,2,3)
imshow(Gauss*par_sol, [min(I(:)), max(I(:))]+3*sig_SNR*[-1,1]) ;
title('Retrieved pattern') ;
subplot(2,2,4)
imshow(C.computeW_(par_sol), []) ;
title('Equivalent weights') ;


