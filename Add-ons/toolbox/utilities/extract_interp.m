%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to extract a region of interest (ROI) in a picture according to
% a set of parameters
% 
% Created: 04/27/2015 (mm/dd/yyyy) Inspired from extract_matching of my PhD
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #pic_in# the picture in which the ROI must be extracted
%
% #[y_c, x_c]# center of the ROI in the picture
%
% #theta# the rotation angle (deg) of the ROI
%
% #[nb_y, nb_x]# size of the wanted ROI
%
% #extrapVal# extrapolation value (default 0)
%
% #scale_factor# factor to scale the image (default 1)
%
% #interp_method# the interpolation method
%   -> 'nearest' - nearest neighbor interpolation
%   -> 'linear'  - bilinear interpolation
%   -> 'spline'  - spline interpolation
%   -> 'cubic'   - bicubic interpolation as long as the data is
%       uniformly spaced, otherwise the same as 'spline' (default)
%
%%%%%%%%
% Ouput
%%%%%%%%
% #pic_out# the extracter ROI
%   /!\ The default value of the extracted pictures is 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pic_out = extract_interp(pic_in, y_c, x_c, theta, nb_y, ...
    nb_x, extrapVal, scale_factor, interp_method)

    %% Initialization
    if nargin < 7
        extrapVal = 0 ;
    end
    if nargin < 8
        scale_factor = 1 ;
    end
    if nargin < 9
        interp_method = 'cubic' ;
    end
    
    pic_in = double(pic_in) ;
    
    % Coordinates of the x-axis
    x_ext = get_Fourier_vector(nb_x, 1) ;
    
    % Coordinates of the y-axis
    y_ext = get_Fourier_vector(nb_y, 1) ;
    
    [x_ext, y_ext] = meshgrid(scale_factor*x_ext, scale_factor*y_ext) ;
    
    %% Interpolation
    if theta~=0
        pic_out = interp2(pic_in, ...
            cosd(theta_rad) * x_ext + sind(theta_rad) * y_ext  + x_c, ...
            -sind(theta_rad) * x_ext + cosd(theta_rad) * y_ext + y_c, ...
            interp_method, extrapVal) ;
    else
        pic_out = interp2(pic_in, x_ext + x_c, y_ext + y_c, ...
            interp_method, extrapVal) ;
    end
end
