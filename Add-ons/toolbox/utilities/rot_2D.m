%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to rotate a pair of matrix.
% 
% Created: 04/10/2018 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #theta# the rotation angle (deg)
%
% #[x_in, y_in]# the input vectors to rotate of the angle theta
%
%%%%%%%%
% Ouput
%%%%%%%%
% #[x_out, y_out]# the output vectors rotated of the angle theta
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_out, y_out] = rot_2D(theta, x_in, y_in)
    x_out = cosd(theta).*x_in - sind(theta).*y_in ;
    y_out = sind(theta).*x_in + cosd(theta).*y_in ;
end