%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to reduce convert a matrix to a gray scale. This matrix can be
% monochromatic or in color.
% 
% Created: 01/14/2015 (mm/dd/yyyy)
% Modified: 12/03/2015 (mm/dd/yyyy) (Suppression of the axial rotation
% parameters)
% Modified: 11/02/2016 (mm/dd/yyyy) (Suppression of the call to mat2gray)
% Modified: 04/25/2017 (mm/dd/yyyy) (Call to rgb2gray and possibility to
%   remove the normalization step)
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #mat_in# the input matrix
%
% #flag_norm# normalization of the input
%   0 -> no normalization (default)
%   1 -> normalization by the mean
%   2 -> normalization by the median
%   3 -> normalization to have the output spreading over the range [0,1]
%
%%%%%%%%
% Ouput
%%%%%%%%
% #mat_out# the output matrix 
%   /!\ The picture is normalized to one
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mat_out = get_gray(mat_in, flag_norm)

    if size(mat_in,3) > 1
        mat_out = double(rgb2gray(mat_in)) ;
    else
        mat_out = double(mat_in) ;
    end
    
    if nargin > 1
        switch flag_norm
            case 1
                mat_out = mat_out./mean(mat_out(:)) ;
            case 2
                mat_out = mat_out./median(mat_out(:)) ;
            case 3
                mx = max(mat_out(:)) ;
                mn = min(mat_out(:)) ;
                mat_out = (mat_out-mn)/(mx-mn) ;
        end
    end
    
end