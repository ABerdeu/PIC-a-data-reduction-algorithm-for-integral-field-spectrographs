% This macro is a procedure to test the LinOpInterp operator
%
% This file supposes that the following toolboxes are in the Matlab path:
%   -> GlobalBioIm
%   -> add-on_GlobalBioIm
%
% Created: 04/05/2018 (mm/dd/yyyy)
% Modified: 09/17/2018 (mm/dd/yyyy) Git desposit
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%


clear all ;
close all ;
clc ;


%% Test 1D
% Parameters
nb_x = 200 ;
nb_q = 75 ; % Number of random query points

% Test
disp('Test 1D') ;
x = 100*rand(nb_x, 1) ;
x = sort(x) ;

x_q = 110*rand(nb_q, 1) ;
x_q = sort(x_q) ;

% x_q = x ;

V = x.*sin(2*pi/25*x) ;

disp('Construction:') ;
disp('LinOpInterp') ;
tic
Interp = LinOpInterp({x}, x_q, 0) ;
disp([num2str(toc), ' sec']) ;
% % disp('LinOpInterp_old') ;
% % tic
% % Interp_old = LinOpInterp_old({x}, x_q, 0) ;
% % disp([num2str(toc), ' sec']) ;

disp('Computation:') ;
disp('LinOpInterp') ;
tic
V_q = Interp*V ;
V_t = Interp'*V_q ;
disp([num2str(toc), ' sec']) ;
% % disp('LinOpInterp_old') ;
% % tic
% % V_q_old = Interp_old*V ;
% % V_t_old = Interp_old'*V_q_old ;
% % disp([num2str(toc), ' sec']) ;

figure(1)
hold on
plot(x,V,'r') ;
plot(x_q,V_q,'b') ;
plot(x,V_t,'g') ;
hold off
legend('theory', 'interp', 'transpos') ;


% % disp('Differences') ;
% % disp(['Interp: ', num2str(mean(V_q-V_q_old))]) ;
% % aux = V_t-V_t_old ;
% % disp(['Interp'': ', num2str(mean(aux(:)))]) ;

checkLinOp(Interp,1) ;
disp([newline, newline]) ;

%% Test 2D
% Parameters
nb_x_2D = 300 ;
nb_y_2D = 450 ;
nb_q_2D = 300 ; % Number of random query points
delta_edge = 0 ; % Extension outward the known values to test the
    % extrapolation technique

% Test
x_2D = 100*rand(nb_x_2D, 1) ;
x_2D = sort(x_2D) ;
y_2D = 75*rand(nb_y_2D, 1) ;
y_2D = sort(y_2D) ;
[xx_2D, yy_2D] = meshgrid(x_2D, y_2D) ;


x_q_2D = (max(x_2D)+2*delta_edge-min(x_2D))*rand(nb_q_2D, 1)+min(x_2D) ...
    - delta_edge ;
y_q_2D = (max(y_2D)+2*delta_edge-min(y_2D))*rand(nb_q_2D, 1)+min(y_2D) ...
    - delta_edge ;

V = 10*xx_2D-3*yy_2D ;
V_th = 10*x_q_2D-3*y_q_2D ;


disp('Test 2D') ;
disp('Construction:') ;
disp('LinOpInterp') ;
tic
Interp = LinOpInterp({y_2D, x_2D}, [y_q_2D, x_q_2D], 0) ;
disp([num2str(toc), ' sec']) ;
% % disp('LinOpInterp_old') ;
% % tic
% % Interp_old = LinOpInterp_old({y_2D, x_2D}, [y_q_2D, x_q_2D], 0) ;
% % disp([num2str(toc), ' sec']) ;

disp('Computation:') ;
disp('LinOpInterp') ;
tic
V_q = Interp*V ;
disp([num2str(toc), ' sec']) ;
tic
V_t = Interp'*V_q ;
disp([num2str(toc), ' sec']) ;
% % disp('LinOpInterp_old') ;
% % tic
% % V_q_old = Interp_old*V ;
% % disp([num2str(toc), ' sec']) ;
% % tic
% % V_t_old = Interp_old'*V_q_old ;
% % disp([num2str(toc), ' sec']) ;


% % disp('Differences Interp - Interp_old') ;
% % disp(['Interp: ', num2str(mean(V_q-V_q_old))]) ;
% % aux = V_t-V_t_old ;
% % disp(['Interp'': ', num2str(mean(aux(:)))]) ;

disp('Differences Interp - Theory') ;
disp(['mean: ', num2str(mean(V_q-V_th))]) ;
disp(['max: ', num2str(max(V_q-V_th))]) ;


figure(2)
hold on
s1 = surf(xx_2D,yy_2D,V) ;
s1.EdgeColor = 'none';
s2 = scatter3(x_q_2D,y_q_2D,V_q,'filled') ;
s3 = surf(xx_2D,yy_2D,V_t) ;
s3.EdgeColor = 'none';
s3.FaceAlpha = 0.5 ;
hold off

checkLinOp(Interp,1) ;
disp([newline, newline]) ;

%% Test 3D
% Parameters
nb_x_3D = 350 ;
nb_y_3D = 400 ;
nb_z_3D = 450 ;
nb_q_3D = 100 ; % Number of random query points

% test
x_3D = 10*rand(nb_x_3D, 1) ;
x_3D = sort(x_3D) ;
y_3D = 5*rand(nb_y_3D, 1) ;
y_3D = sort(y_3D) ;
z_3D = 3*rand(nb_z_3D, 1) ;
z_3D = sort(z_3D) ;
[xx_3D, yy_3D, zz_3D] = meshgrid(x_3D, y_3D, z_3D) ;

x_q_3D = (max(x_3D)-min(x_3D))*rand(nb_q_3D, 1)+min(x_3D) ;
y_q_3D = (max(y_3D)-min(y_3D))*rand(nb_q_3D, 1)+min(y_3D) ;
z_q_3D = (max(z_3D)-min(z_3D))*rand(nb_q_3D, 1)+min(z_3D) ;


V = 10*xx_3D-5*yy_3D+50*zz_3D ;
V_th = 10*x_q_3D-5*y_q_3D+50*z_q_3D ;

disp('Test  3D') ;
disp('Construction:') ;
disp('LinOpInterp') ;
tic
Interp = LinOpInterp({y_3D, x_3D, z_3D}, [y_q_3D, x_q_3D, z_q_3D]) ;
disp([num2str(toc), ' sec']) ;
disp('LinOpInterp_old') ;
% % tic
% % Interp_old = LinOpInterp_old({y_3D, x_3D, z_3D}, [y_q_3D, x_q_3D, z_q_3D]) ;
% % disp([num2str(toc), ' sec']) ;

disp('Computation:') ;
disp('LinOpInterp') ;
tic
V_q = Interp*V ;
disp([num2str(toc), ' sec']) ;
tic
V_t = Interp'*V_q ;
disp([num2str(toc), ' sec']) ;
disp('LinOpInterp_old') ;
% % tic
% % V_q_old = Interp_old*V ;
% % disp([num2str(toc), ' sec']) ;
% % tic
% % V_t_old = Interp_old'*V_q_old ;
% % disp([num2str(toc), ' sec']) ;

% % disp('Differences Interp - Interp_old') ;
% % disp(['Interp: ', num2str(mean(V_q-V_q_old))]) ;
% % aux = V_t-V_t_old ;
% % disp(['Interp'': ', num2str(mean(aux(:)))]) ;

disp('Differences Interp - Theory') ;
disp(['mean: ', num2str(mean(V_q-V_th))]) ;
disp(['max: ', num2str(max(V_q-V_th))]) ;

checkLinOp(Interp,1) ;


%% Test 3D/2D
% Parameters
delta_plot = 5 ; % Shift on the z axis to plot the comparison
ind = 20 ; % Index on which the test is performed
nb_x_3D = 350 ;
nb_y_3D = 400 ;
nb_z_3D = 450 ;
nb_q_3D = 100 ; % Number of random query points


% Test
x_3D = 10*rand(nb_x_3D, 1) ;
x_3D = sort(x_3D) ;
y_3D = 5*rand(nb_y_3D, 1) ;
y_3D = sort(y_3D) ;
z_3D = 3*rand(nb_z_3D, 1) ;
z_3D = sort(z_3D) ;
[xx_3D, yy_3D, zz_3D] = meshgrid(x_3D, y_3D, z_3D) ;

x_q_3D = (max(x_3D)-min(x_3D))*rand(nb_q_3D, 1)+min(x_3D) ;
y_q_3D = (max(y_3D)-min(y_3D))*rand(nb_q_3D, 1)+min(y_3D) ;
z_q_3D = (max(z_3D)-min(z_3D))*rand(nb_q_3D, 1)+min(z_3D) ;

V = 10*xx_3D-5*yy_3D+50*zz_3D ;

disp('Test  XD/3D') ;
disp('Construction:') ;
tic
Interp = LinOpInterp({y_3D, x_3D}, [y_q_3D, x_q_3D], [], ...
    [nb_y_3D, nb_x_3D, nb_z_3D], [1,2]) ;
Interp_2D = LinOpInterp({y_3D, x_3D}, [y_q_3D, x_q_3D]) ;


disp('Computation:') ;
disp('    3D:') ;
tic
V_q = Interp*V ;
V_t = Interp'*V_q ;
disp([num2str(toc), ' sec']) ;

disp('    2D:') ;
tic
V_q_2D = Interp_2D*squeeze(V(:,:,ind)) ;
V_t_2D = Interp_2D'*V_q_2D ;
disp([num2str(toc), ' sec']) ;

X = x_3D;
Y = y_3D ;
X_q = x_q_3D ;
Y_q = y_q_3D ;
[XX_2D, YY_2D] = meshgrid(X, Y) ;

figure(4) ;
hold on
s1 = surf(XX_2D,YY_2D,squeeze(V(:,:,ind))) ;
s1.EdgeColor = 'none';
scatter3(X_q, Y_q, V_q(:, ind),'filled') ;
scatter3(X_q, Y_q, V_q_2D+delta_plot,'filled') ;
s3 = surf(XX_2D,YY_2D,squeeze(V_t(:,:,ind))) ;
s3.EdgeColor = 'none';
s3.FaceAlpha = 0.5 ;
s3 = surf(XX_2D,YY_2D,squeeze(V_t_2D)+5) ;
s3.EdgeColor = 'none';
s3.FaceAlpha = 0.5 ;
hold off

disp('Difference 3D / 2D:') ;
aux = squeeze(V_t(:,:,ind))-squeeze(V_t_2D) ;
disp(max(abs(aux(:)))) ;

% tic
% Interp = LinOpInterp({z_3D}, [z_q_3D], [], ...
%     [nb_y_3D, nb_x_3D, nb_z_3D], [3]) ;
% disp([num2str(toc), ' sec']) ;
% 
% disp('Computation:') ;
% tic
% V_q = Interp*V ;
% disp([num2str(toc), ' sec']) ;
% tic
% V_t = Interp'*V_q ;
% disp([num2str(toc), ' sec']) ;

% X = z_3D ;
% X_q = z_q_3D ;
% 
% figure(5) ;
% ind = 35 ;
% ind_2 = 75 ;
% hold on
% plot(X, squeeze(V(ind,ind_2,:)),'r') ;
% plot(X_q, V_q(:,ind,ind_2)+1,'g') ;
% plot(X, squeeze(V_t(ind,ind_2,:))-1,'b') ;
% hold off

% checkLinOp(Interp,1) ;