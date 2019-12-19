%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get the ROI extraction around a spectrum
% 
% Created: 03/07/2019 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #list_lambda# list of the wavelengths in the model
% 
% #[coef_pol_y, coef_pol_x]# coefficients of the polynomial laws
%
% #rad_ROI# radius of the region of interest to extract
%
% #pix_IFS# pixel characteristics of the sensor
%
%%%%%%%%
% Ouput
%%%%%%%%
% #[list_i, list_j]# List of the pixels to extract
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [list_i, list_j] = ...
    get_ROI(list_lambda, coef_pol_y, coef_pol_x, rad_ROI, pix_IFS)

    %% Initialization
    nb_lambda = length(list_lambda) ;
    deg_pol_y = length(coef_pol_y)-1 ;
    deg_pol_x = length(coef_pol_x)-1 ;
    
    %% Extraction of the spectrum positions
    % y-law
    Pol = zeros(nb_lambda, deg_pol_y+1) ;
    for d = 0:deg_pol_y
        Pol(:,d+1) = list_lambda.^d ;
    end
    y_spec = Pol*coef_pol_y ;


    % Estimation of the polynomial x-law
    Pol = zeros(nb_lambda, deg_pol_x+1) ;
    for d = 0:deg_pol_x
        Pol(:,d+1) = list_lambda.^d ;
    end
    x_spec = Pol*coef_pol_x ;

    %% Extraction of the lenslet region by extending by rad_ROI
    % Upper left corner
    pos_ind = [min(y_spec)-rad_ROI, min(x_spec)-rad_ROI] ;
    ind_min = pos2ind(pos_ind, pix_IFS) ;
    ind_min = max(ind_min, 1) ;

    % Lower right corner
    pos_ind = [max(y_spec)+rad_ROI, max(x_spec)+rad_ROI] ;
    ind_max = pos2ind(pos_ind, pix_IFS) ;
    ind_max = min(ind_max, [pix_IFS.nb_y, pix_IFS.nb_x]) ;

    % Extraction of the region of interest
    list_i = ind_min(1):ind_max(1) ;
    list_j = ind_min(2):ind_max(2) ;
end

