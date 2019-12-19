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
nb_1 = 3 ;
nb_2 = 5 ;

%% Test
SelectorObj = my_LinOpSelectorObj([nb_y, nb_x, nb_1, nb_2], ...
    [1,2], [2,4]) ;

% Check with reals
checkLinOp(SelectorObj,1) ;

