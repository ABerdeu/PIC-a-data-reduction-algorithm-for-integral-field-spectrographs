%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to simulation a calibration pattern
% 
% Created: 05/24/2018 (mm/dd/yyyy)
% Modified: 07/03/2018 (mm/dd/yyyy) Removal the possibility to have several
% frames
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #pix# pixel chracteristic
%   pix.nb_x -> number of the pixel on the x-axis
%   pix.nb_y -> number of the pixel on the y-axis
%
% #[list_y, list_x]# list of the positions of the gaussian
%   patterns in the figure for each wavelength
%
% #gauss_theta# orientation of the gaussian pattern (deg)
%
% #gauss_par# parameters of the gaussian patterns for each wavelength
%   -> sig_x: gaussian extension on the x-axis
%   -> sig_y: gaussian extension on the y-axis
%
% #gauss_amp# amplitude of the gaussian patterns for each gaussian pattern
%
% #gauss_off# offset of the gaussian patterns for each gaussian pattern
%
% #rad_ROI# radius in pixel of the region of interest
%
% #verbose# flag to display the percentage of the simulation (default =
%   false)
%
%%%%%%%%
% Ouput
%%%%%%%%
% #IFS_sim# the simulated field
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IFS_sim = get_IFS_simulation(pix, list_y, list_x, gauss_theta, ...
    gauss_par, gauss_amp, gauss_off, rad_ROI, verbose)

    %% Initialization
    [nb_g, nb_l] = size(list_y) ;
    IFS_sim = zeros(pix.nb_y, pix.nb_x, pix.nb_f) ;
    
    % Maximal indexes in the case that the matrix are not replicated along
    % pix.nb_f
    max_amp = size(gauss_amp, 3) ;
    max_off = size(gauss_off, 1) ;
    
    for f = 1:pix.nb_f
        IFS_sim(:,:,f) = gauss_off(min(f,max_off)) ;
    end
    if nargin < 8
        verbose = false ;
    end
    
    %% Loop on the patterns
    if verbose
        firstprint = fprintf('IFS simulation') ;
        lastprint = fprintf('...') ;
        delta_per = 1 ;
        per_aux = 0 ;
    end
    for g = 1:nb_g
        %% Percentage
        if verbose
            per = g/nb_g*100 ;
            if per>per_aux
                fprintf(repmat('\b', 1, lastprint));
                lastprint = fprintf([': ', num2str(per_aux), ' %%']);
                per_aux = (floor(per/delta_per)+1)*delta_per ;
            end
        end

        %% Loop on the wavelength
        for l = 1:nb_l
            %% Extraction of the pattern of the picture
            ind_y = (-rad_ROI:rad_ROI)+round(list_y(g,l)) + ...
                floor(pix.nb_y/2+1) ;
            ind_x = (-rad_ROI:rad_ROI)+round(list_x(g,l)) + ...
                floor(pix.nb_x/2+1) ;

            % Insuring the extracted frame is in the picture
            ind_y = unique(max(1, min(pix.nb_y, ind_y))) ;
            ind_x = unique(max(1, min(pix.nb_x, ind_x))) ;
            
            %% Simulation of the gaussian pattern
            p_sim = get_gauss_simulation(ind_y', ind_x', ...
                [list_y(g,l), list_x(g,l)] + ...
                [floor(pix.nb_y/2+1), floor(pix.nb_x/2+1)], ...
                gauss_par(l,:), gauss_theta) ;

            %% Update of IFS_sim
            for f = 1:pix.nb_f
                IFS_sim(ind_y, ind_x, f) = IFS_sim(ind_y, ind_x, f) + ...
                    gauss_amp(g, l, min(f,max_amp))*p_sim ;
            end
        end
    end
    
    if verbose
        fprintf(repmat('\b', 1, lastprint)) ;
        fprintf(repmat('\b', 1, firstprint)) ;
        disp(' ') ;
    end
end