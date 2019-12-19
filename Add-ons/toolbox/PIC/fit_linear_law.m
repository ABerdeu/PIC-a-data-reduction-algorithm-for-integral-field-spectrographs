%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to fit a linear law pixelwise.
% 
% Created: 03/03/2019 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #x# the known x on the x-axis, identical for each pixel
%
% #y# the known y on the y-axis on each pixel
%
% #W# Weighting coefficients
%
% #flag_offset# (default: true)
%   true -> y = ax+b
%   false -> y = ax+0
%
%%%%%%%%
% Ouput
%%%%%%%%
% #a# the slop of the fit
%
% #b# the offset of the fit
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a, b] = fit_linear_law(x, y, W, flag_offset)

    %%  Initialization
    if nargin < 3 || isempty(W)
        W = ones(size(y)) ;
    end
    if nargin < 4
        flag_offset = true ;
    end
    
    
    nb_x = length(x) ;
    x = reshape(x, [1,1,nb_x]) ;
    sxy = 0 ;
    for i = 1:nb_x
        sxy = sxy+W(:,:,i).*y(:,:,i)*x(i) ;
    end
    sx = sum(W.*x, 3) ;
    sx2 = sum(W.*x.^2, 3) ;
    W_s = sum(W,3) ;
    
    if flag_offset
        sy = sum(W.*y, 3) ;
        norm = 1./(W_s.*sx2-sx.^2) ;
        a = norm.*(W_s.*sxy-sx.*sy) ;
        b = norm.*(sx2.*sy-sx.*sxy) ;
    else
        a = sxy./sx2 ;
        b = 0 * a ;
    end
end