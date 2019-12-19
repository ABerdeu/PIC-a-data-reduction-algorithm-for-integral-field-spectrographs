% This macro is a procedure to test the my_LinOpSelectorPatch operator
%
% This file supposes that the following toolboxes are in the Matlab path:
%   -> GlobalBioIm
%   -> add-on_GlobalBioIm
%
% Created: 06/05/2018 (mm/dd/yyyy)
% Modified: 09/17/2018 (mm/dd/yyyy) Git desposit
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%


clear all ;
close all ;
clc ;


%% Parameters
nb_x = 200 ;
nb_y = 450 ;
nb_z = 100 ;
flag_squeeze = true ;

%% Test
SelectorPatch = my_LinOpSelectorPatch([nb_y, nb_x, nb_z], ...
    [2,5,8], [2,10,8], flag_squeeze) ;

disp(SelectorPatch.sizeout) ;
disp(SelectorPatch*rand([nb_y, nb_x, nb_z])) ;


% Check with reals
checkLinOp(SelectorPatch,1) ;

