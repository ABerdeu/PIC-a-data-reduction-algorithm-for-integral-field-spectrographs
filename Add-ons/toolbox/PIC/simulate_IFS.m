%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to simulate an IFS for a given set of calibration parameters for
% symmetric patterns
% 
% Created: 09/04/2018 (mm/dd/yyyy)
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
% #pattern_model# flags to specify the pattern to fit
%   pattern_model.flag_norm -> Normalized?
%   pattern_model.flag_profile -> Gaussian / Moffat
%   pattern_model.oversampling -> oversampling of the pattern (default: 1)
%
% #[list_y, list_x]# list of the positions of the patterns in the figure
%   for each wavelength
%
% #par_pat# parameters of the pattern
% 
% #list_amp# of the amplitude of the pattern's realization in the picture
%
% #offset# offset of the patterns
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
function IFS_sim = simulate_IFS(pix, pattern_model, list_y, ...
    list_x, par_pat, list_amp, offset, rad_ROI, verbose)

    %% Initialization
    [nb_pat, nb_l] = size(list_y) ;
    IFS_sim = zeros(pix.nb_y, pix.nb_x, pix.nb_f) ;
    
    % Maximal indexes in the case that the matrix are not replicated along
    % pix.nb_f
    max_amp = size(list_amp, 3) ;
    max_off = size(offset, 1) ;
    for f = 1:pix.nb_f
        IFS_sim(:,:,f) = offset(min(f,max_off)) ;
    end
    if nargin < 8
        verbose = false ;
    end
    
    % Maximal indexes for the pattern parameters (are they the same for all
    % patterns or different?)
    if size(par_pat, 1)==1
        max_par = 1 ;
    else
        max_par = nb_pat ;
    end
    
    %% Loop on the patterns
    if verbose
        delta_per = 1 ;
        [per_aux, lastprint] = display_percentage('init', ...
            'IFS simulation') ;
    end
    for pat = 1:nb_pat
        %% Percentage
        if verbose
            [per_aux, lastprint] = display_percentage('iter', ...
                {pat/nb_pat*100, per_aux, delta_per, lastprint}) ;
        end

        %% Loop on the wavelength
        for l = 1:nb_l
            %% Extraction of the pattern of the picture
            ind_y = (-rad_ROI:rad_ROI)+round(list_y(pat,l)) + ...
                floor(pix.nb_y/2+1) ;
            ind_x = (-rad_ROI:rad_ROI)+round(list_x(pat,l)) + ...
                floor(pix.nb_x/2+1) ;

            % Insuring the extracted frame is in the picture
            ind_y = unique(max(1, min(pix.nb_y, ind_y))) ;
            ind_x = unique(max(1, min(pix.nb_x, ind_x))) ;
            
            %% Simulation of the gaussian pattern
            par_pat_l = get_pattern_param(par_pat, l, min(pat,max_par)) ;

            % Calibration pattern in the image corner
            p_sim = get_pattern_simulation(ind_y', ind_x', ...
                list_y(pat,l) + floor(pix.nb_y/2+1), ...
                list_x(pat,l) + floor(pix.nb_x/2+1), par_pat_l, ...
                1, 0, pattern_model) ;

            %% Update of IFS_sim
            for f = 1:pix.nb_f
                IFS_sim(ind_y, ind_x, f) = IFS_sim(ind_y, ind_x, f) + ...
                    list_amp(pat, l, min(f,max_amp))*p_sim ;
            end
        end
    end
    display_percentage('exit', lastprint) ;
end