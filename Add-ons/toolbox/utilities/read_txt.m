%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to read a 2D matrix in a txt file
% 
% Created: 08/30/2019 (mm/dd/yyyy)
% Author: Anthony Berdeu (Centre de Recherche en Astrophysique de Lyon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #var# the variable to save
%
% #name# the name of the projected file
%
% #path_save# the path to the results
%
%%%%%%%%
% Ouput
%%%%%%%%
% #var# the opend variable
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [var] = read_txt(name, path_save)
    var = load([path_save, name, '.txt'], '-ascii') ;
end