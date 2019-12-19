%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get the coordinates in the Fourier space.
% Because of the fft and the fftshift properties, a zero must be in the
% list at a specific position.
% 
% Created: 03/19/2015 (mm/dd/yyyy)
% Modified: 10/26/2015 (mm/dd/yyyy) (changing the name, from
% "get_kernel_vector")
% Author: Anthony Berdeu (PhD student at CEA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #nb_el# the number of elements in the vector
%
% #d_el# the elementary length of the elements in the vector
%
%%%%%%%%
% Ouput
%%%%%%%%
% #pos_Fourier# the position of the kernel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pos_Fourier] = get_Fourier_vector(nb_el, d_el)

    if mod(nb_el,2)
        % There is an odd number of elements, the center of the frame is at
        % the center of the vector
        pos_Fourier = ((0:nb_el-1)-(nb_el-1)/2)*d_el ;
    else
        % There is an even number of elements, the center of the frame is
        % at the position round(nb_pix/2)+1
        pos_Fourier = ((0:nb_el-1)-nb_el/2)*d_el ;
    end
end