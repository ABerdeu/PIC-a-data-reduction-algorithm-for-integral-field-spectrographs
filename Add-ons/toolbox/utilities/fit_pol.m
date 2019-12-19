%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to fit a polynomial expression in a could of points.
% 
% Created: 08/08/2018 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #deg_pol# degree of the polynomial to fit
%
% #x# the known positions on the x-axis
%
% #y# the known positions on the y-axis on each cloud. It has a size of 
% [number of points of each cloud, number of clouds to fit].
%
% #x_fit# the positions to fit on the x-axis
%
% #[Mat_x, Mat_x_fit]# (optional) the precomputed powers of x and x_fit
%
%%%%%%%%
% Ouput
%%%%%%%%
% #y_fit# the fitted curves [number of points to fit, number of
% clouds to fit]
%
% #coef_pol# coefficient of the polynomial expression
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y_fit, coef_pol] = fit_pol(deg_pol, x, y, x_fit, Mat_x, ...
    Mat_x_fit)

    % Estimating coefficients
    if nargin < 5
        Mat_x = [] ;
    end
    if isempty(Mat_x)
        [Mat_deg, Mat_x] = meshgrid(0:deg_pol, x) ;
        Mat_x = Mat_x.^Mat_deg ;
    end
    coef_pol = Mat_x\y ;

    % Fitting polynomial expression
    if nargin < 6
        Mat_x_fit = [] ;
    end
    if isempty(Mat_x_fit)
        [Mat_deg, Mat_x_fit] = meshgrid(0:deg_pol, x_fit) ; 
        Mat_x_fit = Mat_x_fit.^Mat_deg ;
    end
    y_fit = Mat_x_fit*coef_pol ;
end