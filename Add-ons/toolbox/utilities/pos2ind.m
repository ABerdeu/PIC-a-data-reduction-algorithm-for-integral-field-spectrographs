%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to convert the a list of positions on a picture to the indexes
% in the matrix
% 
% Created: 05/18/2018 (mm/dd/yyyy)
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
% Modified: 07/15/2018 (mm/dd/yyyy)
% Author: Anthony Berdeu (Centre de Recherche en Astrophysique de Lyon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #list_pos = [list_y, list_x]# list of the positions
%
% #pix# pixel chracteristic
%   pix.nb_x -> number of the pixel on the x-axis
%   pix.nb_y -> number of the pixel on the y-axis
%   pix.dx -> size of the pixel on the x-axis
%   pix.dy -> size of the pixel on the y-axis
%
%%%%%%%%
% Ouput
%%%%%%%%
% #list_ind = [list_i, list_j]# list of the indexes (integer)
%
% #list_ind_float = [list_i, list_j]# list of the indexes (float)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [list_ind, list_ind_float] = pos2ind(list_pos, pix)

    list_ind_float = [list_pos(:,1)/pix.dy+floor(pix.nb_y/2)+1 , ...
        list_pos(:,2)/pix.dx+floor(pix.nb_x/2)+1] ;

    list_ind = round(list_ind_float) ;
end