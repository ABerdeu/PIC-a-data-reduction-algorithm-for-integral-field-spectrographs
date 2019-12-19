%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to remove a directory if it is possible
% 
% Created: 04/03/2017 (mm/dd/yyyy)
% Author: Anthony Berdeu (PhD student at CEA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #path_dir# the directory to remove
%
%%%%%%%%
% Ouput
%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function remove_dir(path_dir)
    try
        rmdir(path_dir, 's') ;
    catch
        disp('Unable to remove directory...') ;
        disp('Directory not removed, manual removal needed...') ;
    end
end