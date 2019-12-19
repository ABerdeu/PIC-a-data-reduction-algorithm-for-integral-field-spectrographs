% This macro is a procedure to test the OpMoffat operator
%
% This file supposes that the following toolboxes are in the Matlab path:
%   -> GlobalBioIm
%   -> add-on_GlobalBioIm
%
% Created: 08/31/2018 (mm/dd/yyyy)
% Modified: 09/17/2018 (mm/dd/yyyy) Git desposit
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%


clear all ;
close all ;
clc ;


%% Parameters
x = (-15:15)' ;
y = (-10:10)' ;
y_c = 3 ;
x_c = -5 ;
alpha = 2 ;
beta = 3 ;
amp = 1 ;
amp_norm = amp/(beta-1)*pi*alpha^2 ; % If normalized...
offset = 0.25 ;


%% On position, amplitude and offset
delta_0 = [1, 0.5, 2, -2]' ; % differences for the initialization of the
% optimization
sig_SNR = 0.1 ;
opti_method = 'Huber' ; % 'Cauchy', 'Huber', ...

% Simulation
Moffat = OpMoffat(y, x, [], [], alpha, beta, [], [], false) ;

I = Moffat*[y_c; x_c; amp; offset] ;
I_s = I+sig_SNR*randn(size(I)) ;
I_s=imnoise(I_s/(amp+offset), 'salt & pepper', 0.1) ;
I_s=I_s*(amp+offset) ;


checkMap(Moffat) ;

% Test optimization
s = max(1.4826*median(abs(I_s(:)-median(I_s(:)))),0.002) ;
% s = 1 ;
options.method = opti_method ;
options.flag_s = s ;
C = CostRobustPenalization(Moffat, I_s, options) ;
const_min = [-Inf; -Inf; 0.01; -Inf] ;
const_max = [+Inf; +Inf; +Inf; +Inf] ;
VMLMB=OptiVMLMB(C,const_min,const_max);  
VMLMB.ItUpOut=1; 
VMLMB.verbose=false ;
VMLMB.maxiter=20000;      % max number of iterations
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
imshow(Moffat*par_sol, []) ;
title('Retrieved pattern') ;
subplot(2,2,4)
imshow(C.computeW_(par_sol), []) ;
title('Equivalent weights') ;

%% On alpha, beta, amp, offset
delta_0 = [1, 0.5, 2, -2]' ; % differences for the initialization of the
% optimization
sig_SNR = 0.01 ;
opti_method = 'Cauchy' ; % 'Cauchy', 'Huber', ...

% Simulation
Moffat = OpMoffat(y, x, y_c, x_c, [], [], [], [], true) ;

I = Moffat*[alpha; beta; amp_norm; offset] ;
I_s = I+sig_SNR*randn(size(I)) ;
I_s=imnoise(I_s/(amp+offset), 'salt & pepper', 0.1) ;
I_s=I_s*(amp+offset) ;

% Test optimization
s = 1.4826*median(abs(I_s(:)-median(I_s(:)))) ;
options.method = opti_method ;
options.flag_s = s ;
C = CostRobustPenalization(Moffat, I_s, options) ;
const_min = [0.1; 0.1; 0; -Inf] ;
const_max = [+Inf; +Inf; +Inf; +Inf] ;
VMLMB=OptiVMLMB(C,const_min,const_max);  
VMLMB.ItUpOut=1; 
VMLMB.maxiter=100000;      % max number of iterations
VMLMB.m=3;              % number of memorized step in hessian approximation
par_th = [alpha; beta; amp_norm; offset] ;
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
imshow(Moffat*par_sol, []) ;
title('Retrieved pattern') ;
subplot(2,2,4)
imshow(C.computeW_(par_sol), []) ;
title('Equivalent weights') ;

%% On position, amplitude, aplha and beta
delta_0 = [0.75, -0.5, 1, 1, 5]' ; % differences for the initialization of
% the optimization
sig_SNR = 0.01 ;
opti_method = 'Cauchy' ; % 'Cauchy', 'Huber', ...

% Simulation
Moffat = OpMoffat(y, x, [], [], [], [], [], offset, true) ;

I = Moffat*[y_c; x_c; alpha; beta; amp_norm] ;
I_s = I+sig_SNR*randn(size(I)) ;
I_s=imnoise(I_s/(amp+offset), 'salt & pepper', 0.1) ;
I_s=I_s*(amp+offset) ;

checkMap(Moffat) ;

% Test optimization
s = max(1.4826*median(abs(I_s(:)-median(I_s(:)))),0.002) ;
% s = 1 ;
options.method = opti_method ;
options.flag_s = s ;
C = CostRobustPenalization(Moffat, I_s, options) ;
const_min = [-Inf; -Inf; 0.1; 0.1; 0] ;
const_max = [+Inf; +Inf; +Inf; +Inf; +Inf] ;
VMLMB=OptiVMLMB(C,const_min,const_max);  
VMLMB.ItUpOut=1; 
VMLMB.verbose=false ;
% VMLMB.maxiter=5;      % max number of iterations
VMLMB.m=3;              % number of memorized step in hessian approximation
par_th = [y_c; x_c; alpha; beta; amp_norm] ;
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
imshow(Moffat*par_sol, []) ;
title('Retrieved pattern') ;
subplot(2,2,4)
imshow(C.computeW_(par_sol), []) ;
title('Equivalent weights') ;
