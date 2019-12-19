%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to save a matrix to a fits file
% 
% Created: 03/20/2015 (mm/dd/yyyy)
% Author: Anthony Berdeu (PhD student at CEA)
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
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = save_fits(var, name, path_save)
    if ~exist(path_save,'dir')
        mkdir(path_save)
    end
    
    path_tot = [path_save, name, '.fits'] ;
    
    if exist(path_tot,'file')
        delete(path_tot) ;
    end
    fitswrite(var,path_tot);
end