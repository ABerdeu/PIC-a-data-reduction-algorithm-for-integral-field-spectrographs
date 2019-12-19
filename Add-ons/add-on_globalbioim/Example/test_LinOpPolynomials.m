% This macro is a procedure to test the LinOpPolynomial operator
%
% This file supposes that the following toolboxes are in the Matlab path:
%   -> GlobalBioIm
%   -> add-on_GlobalBioIm
%
% Created: 10/08/2018 (mm/dd/yyyy)
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%


clear all ;
close all ;
clc ;


%% Test 1D
disp('Test 1D') ;

% Parameters
nb_x = 30 ;
sig = 1 ;
option_opti.method = 'fminunc' ;
option_opti.verbose = false ;
option_opti.maxiter = 1000 ;
degree = 3 ;

% Test
x = ((1:nb_x)-floor(nb_x/2+1))' ;
v = 10 - 5.*x + 0.5.*x.^2 + 0.025.*x.^3 ;
v_n = v + sig^2*randn(size(v)) ;

Pol = LinOpPolynomial(x, degree) ;


C = CostL2(Pol.sizeout,v_n)*Pol ;
[par_out] = run_Opti(C, [], [], option_opti, zeros(Pol.sizein)) ;
v_fit = Pol*par_out ;

disp('Found coefficients:') ;
disp(par_out) ;

fig = figure(1) ;
clf(fig) ;
hold on
plot(x,v,'r') ;
plot(x,v_n,'b') ;
plot(x,v_fit,'g') ;
hold off
legend('theory', 'noisy', 'found') ;


checkLinOp(Pol,1) ;
disp([newline, newline]) ;

%% Test 2D
disp('Test 2D') ;

% Parameters
nb_x = 30 ;
nb_y = 50 ;
sig = 5 ;
option_opti.method = 'VMLMB' ;
option_opti.verbose = false ;
option_opti.maxiter = 1000 ;
degree = 2 ;
delta_plot = 100 ;

% Test
x = ((1:nb_x)-floor(nb_x/2+1))' ;
y = ((1:nb_y)-floor(nb_y/2+1))' ;
[x,y] = meshgrid(x,y) ;
v = 10 - 5.*x + 3.*y + 0.5.*x.^2 - 0.25.*y.^2 + 0.75.*y.*x ;
v_n = v + sig^2*randn(size(v)) ;

Pol = LinOpPolynomial({y,x}, degree) ;


C = CostL2(Pol.sizeout,v_n)*Pol ;
[par_out] = run_Opti(C, [], [], option_opti, zeros(Pol.sizein)) ;
v_fit = Pol*par_out ;

disp('Found coefficients:') ;
disp(par_out) ;

fig = figure(2) ;
clf(fig) ;
hold on
s1 = surf(x,y,v) ;
s1.EdgeColor = 'none';
s2 = surf(x,y,v_n+delta_plot) ;
s2.EdgeColor = 'none';
s2.FaceAlpha = 0.5 ;
s3 = surf(x,y,v_fit+2*delta_plot) ;
s3.EdgeColor = 'none';
s3.FaceAlpha = 0.5 ;
hold off

checkLinOp(Pol,1) ;
disp([newline, newline]) ;

%% Test 3D
disp('Test 3D') ;

% Parameters
nb_x = 30 ;
nb_y = 50 ;
nb_z = 25 ;
sig = 10 ;
option_opti.method = 'VMLMB' ;
option_opti.verbose = false ;
option_opti.maxiter = 1000 ;
degree = 2 ;
delta_plot = 100 ;
z_plot = round(nb_z/2) ;

% Test
x = ((1:nb_x)-floor(nb_x/2+1))' ;
y = ((1:nb_y)-floor(nb_y/2+1))' ;
z = ((1:nb_z)-floor(nb_z/2+1))' ;
[x,y,z] = meshgrid(x,y,z) ;
v = 10 - 5.*x + 3.*y + 0.5.*x.^2 - 0.25.*y.^2 + 0.75.*y.*x + ...
    0.85.*x.*z  - 0.65.*y.*z - 0.35.*z.^2 - 3.*z ;
v_n = v + sig^2*randn(size(v)) ;

Pol = LinOpPolynomial({y,x,z}, degree) ;


C = CostL2(Pol.sizeout,v_n)*Pol ;
[par_out] = run_Opti(C, [], [], option_opti, zeros(Pol.sizein)) ;
v_fit = Pol*par_out ;

disp('Found coefficients:') ;
disp(par_out) ;

fig = figure(3) ;
clf(fig) ;
hold on
s1 = surf(x(:,:,z_plot),y(:,:,z_plot),v(:,:,z_plot)) ;
s1.EdgeColor = 'none';
s2 = surf(x(:,:,z_plot),y(:,:,z_plot),v_n(:,:,z_plot)+delta_plot) ;
s2.EdgeColor = 'none';
s2.FaceAlpha = 0.5 ;
s3 = surf(x(:,:,z_plot),y(:,:,z_plot),v_fit(:,:,z_plot)+2*delta_plot) ;
s3.EdgeColor = 'none';
s3.FaceAlpha = 0.5 ;
hold off

checkLinOp(Pol,1) ;
disp([newline, newline]) ;