%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to fit the positions of a spectrum in calibration points.
% 
% Created: 08/08/2018 (mm/dd/yyyy)
% Modified: 10/18/2018 (mm/dd/yyyy) Better linear fit
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #deg_pol_spec=[deg_pol_spec_xl, deg_pol_spec_y]# Degree
%   of the polynomial to fit the spectra position along the dispersion
%   direction(xl) and its orthogonal direction (y) obtained with the
%   calibration data
%
% #pos_cal = [x_cal, y_cal]# the calibration points
%
% #list_lambda_cal# the list of the calibration wavelengths
%
% #list_lambda# the list of the wavelengths to fit
%
%%%%%%%%
% Ouput
%%%%%%%%
% #[x_fit, y_fit]# the fitted spectrum in the original frame
%
% #[x_fit_local, y_fit_local]# the fitted spectrum in the frame of the
% spectrum
%
% #lin_fit = [origin, theta]# parameters of the linear fit to obtain the
% local frame of the spectrum
%   origin -> origin of the linear fit
%   theta -> orientation of the linear fit (deg)
%
% #[x_cal, y_cal]# positions of the calibration points in the local frame
% of the spectrum
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_fit, y_fit, x_fit_local, y_fit_local, lin_fit, x_cal, ...
    y_cal] = fit_spectrum(deg_pol_spec, pos_cal, list_lambda_cal, ...
list_lambda)

    % Initialization
    x_cal = pos_cal(:,1) ;
    y_cal = pos_cal(:,2) ;
    theta = mod(atand((y_cal(end)-y_cal(1))/(x_cal(end)-x_cal(1))), 180) ;
    theta = mod(theta, 180) ;
    origin = [mean(x_cal), mean(y_cal)] ;

    % Linear fit to find the local frame of the spectrum:
    % lin_fit = [origin = [x_c, y_c] , theta]
    % The distance is given by the cross product
    dist2 = @(lin_fit)sum( ...
        ( ...
        cosd(lin_fit(3)).*(y_cal-lin_fit(2)) ...
        - ...
        sind(lin_fit(3)).*(x_cal-lin_fit(1)) ...
        ).^2) ;
    options.MaxIter = 25000 ;
    lin_fit = run_OptiFMINSEARCH(dist2, [], [], options, ...
        [origin, theta]) ;
    x_cal = x_cal-lin_fit(1) ;
    y_cal = y_cal-lin_fit(2) ;
    
    % Orginal to local frame
    theta_spec = lin_fit(3) ;
    [x_cal, y_cal] = rot_2D(-theta_spec, x_cal, y_cal) ;
    
    % Fitting polynomial expression
    x_fit_local = fit_pol(deg_pol_spec(1), list_lambda_cal, x_cal, ...
        list_lambda) ;
    y_fit_local = fit_pol(deg_pol_spec(2), x_cal, y_cal, x_fit_local) ;
    
    % Going back to original frame
    [x_fit, y_fit] = rot_2D(theta_spec, x_fit_local, y_fit_local) ;
    x_fit = x_fit+lin_fit(1) ;
    y_fit = y_fit+lin_fit(2) ;
end