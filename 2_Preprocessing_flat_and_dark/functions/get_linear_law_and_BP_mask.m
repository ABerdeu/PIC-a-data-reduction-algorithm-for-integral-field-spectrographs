%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get the linear law of a temporal acquisition and the
% resulting bad pixels
% 
% Created: 01/26/2019 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #acqui_list# cell with the list of the acquisitions
%
% #exposure_list# list of the exposures
%
% #tol_MAD# Tolerance in residues for bad pixels indentification
%
%%%%%%%%
% Ouput
%%%%%%%%
% #linear_law# the coefficients of the linear law
%
% #BP_mask# the mask of the bad pixels
%
% #residues# the residues from which the bad pixels are extracted
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [linear_law, BP_mask, residues] = ...
    get_linear_law_and_BP_mask(acqui_list, exposure_list, tol_MAD)

    %% Initialization
    nb_y = size(acqui_list{1}, 1) ;
    nb_x = size(acqui_list{1}, 2) ;
    nb_acqui = length(acqui_list) ;
    AVG_acqui_list = zeros(nb_y, nb_x, nb_acqui) ; % List of the averaged
        % acquisitions
    linear_law = zeros(nb_y, nb_x) ; % Temporaly average acquisition
    for acqui = 1:nb_acqui
        disp(['    Acqui: ', num2str(acqui), '/', num2str(nb_acqui)]) ;
        
        % Extracting averaged acqui
        AVG_acqui_list(:,:,acqui) = mean(acqui_list{acqui}, 3) ;
        
        % Getting linear law by last mean square fit on each pixels
        linear_law = linear_law + ...
            AVG_acqui_list(:,:,acqui)*exposure_list(acqui) ;
    end
    linear_law = linear_law / sum(exposure_list.^2) ;

    %% Analyzing dark residues
    BP_mask = ones(nb_y, nb_x) ;
    residues = zeros(nb_y, nb_x, nb_acqui) ;
    for acqui = 1:nb_acqui
        disp(['    Acqui residues: ', num2str(acqui), '/', ...
            num2str(nb_acqui)]) ;
        residues(:,:,acqui) = AVG_acqui_list(:,:,acqui) - ...
            exposure_list(acqui)*linear_law ;
        
        % MAD Statistics
        stat = residues(:,:,acqui) ;
        med_stat = median(stat(:)) ;
        MAD_stat = median(abs(stat(:)-med_stat)) ;

        % Masking non linear pixels
        BP_mask = BP_mask.*(abs(stat-med_stat)<tol_MAD*MAD_stat) ;
    end
end