%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to open a stack in a tif format
% 
% Created: 07/07/2017 (mm/dd/yyyy)
% Author: Anthony Berdeu (PhD student at CEA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #path_data# path of the tif
%
% #name_data# name of the tif
%
% #flip_y# flipping the y-axis
%   1 -> yes
%   0 -> no
%
% #flip_z# flipping the z-axis
%   1 -> yes
%   0 -> no
%
%%%%%%%%
% Ouput
%%%%%%%%
% #stack# the open stack
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stack = open_tif(path_data, name_data, flip_y, flip_z)

    %% Parameters for the volume extraction
    info_tif = imfinfo([path_data, name_data]) ;
    nb_y = info_tif.Height ;
    nb_x = info_tif.Width ;
    nb_z = length(info_tif) ;
    BitDepth = info_tif.BitDepth ;
    switch BitDepth
        case 8
            type_dat = 'uint8' ;
        case 16
            type_dat = 'uint16' ;
        otherwise
            type_dat = 'double' ;
    end


    %% Loading the volume slice by slice
    if nb_z>1
        stack = zeros(nb_y, nb_x, nb_z, type_dat) ;
        for k = 1:nb_z
            disp(['        Opening slice ', num2str(k), '/', ...
                num2str(nb_z)]) ;
            currentImage = imread([path_data, name_data], k, ...
                'Info', info_tif);
            if ~flip_y % By default, the tif are flipped in matlab
                currentImage = flip(currentImage,1) ;
            end

            if ~flip_z % By default, the tif are flipped in matlab
                stack(:,:,k) = currentImage;
            else
                stack(:,:,nb_z+1-k) = currentImage;
            end
        end
    else
        stack = imread([path_data, name_data], 1, 'Info', info_tif);
    end
end