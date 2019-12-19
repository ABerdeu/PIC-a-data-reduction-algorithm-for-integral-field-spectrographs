%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get the estimation of the variance law in function of the
% median value of adu
% 
% Created: 03/04/2019 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #list_avg# list of the averaged files in adu
%
% #list_var# list of the variance of the files in adu
%
%%%%%%%%
% Ouput
%%%%%%%%
% #list_med# list of the median values on the acquisitions
%
% #list_var_med# list of the estimated variance
%
% #list_var_MAD# list of the robust estimation of the standard deviation
% on the variance
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [list_med, list_var_med, list_var_MAD] = ...
    get_variance_estimate(list_avg, list_var)

    % Loop on the acquisitions
    nb_acqui = size(list_avg, 3) ;
    list_med = zeros(nb_acqui, 1) ;
    list_var_med = zeros(nb_acqui, 1) ;
    list_var_MAD = zeros(nb_acqui, 1) ;
    for acqui = 1:nb_acqui
        disp(['    Acqui: ', num2str(acqui), '/', num2str(nb_acqui)]) ;

        % Median value of the dark
        aux = list_avg(:,:,acqui) ;
        list_med(acqui) = median(aux(:)) ;

        % Analyzing variance
        aux = list_var(:,:,acqui) ;
        list_var_med(acqui) = median(aux(:)) ;
        list_var_MAD(acqui) = ...
            median(abs(aux(:)-list_var_med(acqui))) ;
    end
end