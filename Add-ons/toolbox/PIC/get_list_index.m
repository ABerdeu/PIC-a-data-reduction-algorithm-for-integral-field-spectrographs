%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get a list of index in a picture around a given position for
% a given size
% 
% Created: 07/26/2018 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #pos = [y_c, x_c]# the position of the region of interest
%
% #rad_ROI = [rad_ROI_y, rad_ROI_x]# the radius of the region of interest
%
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
% #[ind_y, ind_x]# the list of the indexes in the matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ind_y, ind_x] = get_list_index(pos, rad_ROI, pix)

    % Local extraction 
    ind_y = (-rad_ROI(1):rad_ROI(1))+round(pos(1)) + ...
        floor(pix.nb_y/2+1) ;
    ind_x = (-rad_ROI(2):rad_ROI(2))+round(pos(2)) + ...
        floor(pix.nb_x/2+1) ;

    % Insuring the extracted frame is in the picture
    ind_y = unique(max(1, min(pix.nb_y, ind_y)))' ;
    ind_x = unique(max(1, min(pix.nb_x, ind_x)))' ;
end