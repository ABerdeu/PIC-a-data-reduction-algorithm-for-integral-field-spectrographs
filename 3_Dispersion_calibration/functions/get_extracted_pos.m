%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get the extracted position in a grid from the given extracted
% indexes
% 
% Created: 03/07/2019 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #[list_i, list_j]# List of the extracted pixels
% 
% #[coef_pol_y, coef_pol_x]# coefficients of the polynomial laws
%
% #pix_IFS# pixel characteristics of the sensor
%
%%%%%%%%
% Ouput
%%%%%%%%
% #[crop_y, crop_x]# List of the extracted positions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [crop_y, crop_x] = get_extracted_pos(list_i, list_j, pix_IFS)

    %% Shift in the positions to change to the new frame
    delta_pos_center = [ ...
        floor(pix_IFS.nb_y/2)+1 - ...
            list_i(floor(length(list_i)/2)+1) ; ...
        floor(pix_IFS.nb_x/2)+1 - ...
            list_j(floor(length(list_j)/2)+1) ] ;
    pix_crop = pix_IFS ;
    pix_crop.nb_y = length(list_i) ;
    pix_crop.nb_x = length(list_j) ;

    %% Positions of the pixels on the sensor
    crop_x = get_Fourier_vector(pix_crop.nb_x, pix_crop.dx) - ... 
        delta_pos_center(2) ;
    crop_y = get_Fourier_vector(pix_crop.nb_y, pix_crop.dy) - ... 
        delta_pos_center(1) ;
end

