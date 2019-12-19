%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to save a 2D matrix to a txt file
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
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = save_txt(var, name, path_save)
    
    if length(size(var)) > 2
        error('The input variable has more than three dimensions...') ;
    else
        if ~exist(path_save,'dir')
            mkdir(path_save)
        end

        path_tot = [path_save, name, '.txt'] ;

        if exist(path_tot,'file')
            delete(path_tot) ;
        end
        save(path_tot, 'var', '-ascii', '-double') ;
    end
end