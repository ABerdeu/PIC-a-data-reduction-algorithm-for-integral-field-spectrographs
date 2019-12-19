% This macro is a procedure to test the LinOpPad operator
%
% This file supposes that the following toolboxes are in the Matlab path:
%   -> GlobalBioIm
%   -> add-on_GlobalBioIm
%
% Created: 05/02/2018 (mm/dd/yyyy)
% Modified: 09/17/2018 (mm/dd/yyyy) Git desposit
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%


clear all ;
close all ;
clc ;


%% Parameters
nb_x = 100 ;
nb_y = 200 ;
% nb_x_pad = 100 ;
% nb_y_pad = 200 ;
nb_x_pad = 120 ;
nb_y_pad = 250 ;

%% Test
% Pad = LinOpPad([nb_y, nb_x], [nb_y_pad, nb_x_pad], 'mean') ; % Padding
%     % the the average value
Pad = LinOpPad([nb_y, nb_x], [nb_y_pad, nb_x_pad]) ; % Padding with the
    % default value

checkLinOp(Pad,1) ;

pic = rand(nb_y, nb_x) ;
pic_pad = Pad*pic ;
pic_2 = Pad'*pic_pad ;

figure(1) ;
subplot(1,3,1) ;
imshow(pic, []) ;
title('pic') ;

subplot(1,3,2) ;
imshow(pic_pad, []) ;
title('Padded pic') ;

subplot(1,3,3) ;
imshow(pic_2, []) ;
title('Trans padded pic') ;