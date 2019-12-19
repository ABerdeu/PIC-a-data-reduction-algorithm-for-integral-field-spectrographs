%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get the projection matrix to project the spectra diffracted
% by the lenslet array on the sensor
% 
% Created: 11/19/2018 (mm/dd/yyyy)
% Created: 11/27/2018 (mm/dd/yyyy) Adaptation to polynomial laws
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #pix# pixel chracteristic
%   pix.nb_x -> number of the pixel on the x-axis
%   pix.nb_y -> number of the pixel on the y-axis
%   pix.dx -> size of the pixel on the x-axis
%   pix.dy -> size of the pixel on the y-axis
%
% #rad_ROI# radius of the pattern (in pixel)
%
% #pattern_model# flags to specify the pattern to fit
%	pattern_model.flag_norm -> Normalized?
%   pattern_model.flag_profile -> Gaussian / Moffat
%   pattern_model.oversampling -> oversampling of the pattern 
%   (default: 1)
%
% #list_lambda# list of the spectra wavelengths
%
% #[list_y_0, list_x_0]# Positions of the absolute grid
%
% #coef_pol_x# list of the coefficients of the polynomial dispersion of the
% spectra on the x-axis
%
% #coef_pol_y# list of the coefficients of the polynomial dispersion of the
% spectra on the y-axis
%
% #coef_pol_dif# list of the coefficients of the polynomial dispersion
%   of the patterns size along the spectra
%
% #deg_pol_x# degrees of the polynomial dispersion of the
% spectra on the x-axis
%
% #deg_pol_y# degrees of the polynomial dispersion of the
% spectra on the y-axis
%
% #deg_pol_dif# degrees of the polynomial dispersion
%   of the patterns size along the spectra
%
% #max_deg_pol_x# maximal degrees of the polynomial dispersion of the
% spectra on the x-axis
%
% #max_deg_pol_y# maximal degrees of the polynomial dispersion of the
% spectra on the y-axis
%
% #max_deg_pol_dif# maximal degrees of the polynomial dispersion
%   of the patterns size along the spectra
%
% #trans# List of the transmission of each lenslet and possibly at all 
%   wavelength (default = 1)
%
% #[correction_y, correction_x] = List of the correction to apply on the
% coefficients of the position law
%
% #correction_dif# List of the correction to apply on the coefficients of
% the diffraction law (default: 0)
%
%%%%%%%%
% Ouput
%%%%%%%%
% #s# the scaling factor
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Spec_proj = get_SpecProjMat(pix, rad_ROI, pattern_model, ...
    list_lambda, list_y_0, list_x_0, coef_pol_x, coef_pol_y, ...
    coef_pol_dif, deg_pol_x, deg_pol_y, deg_pol_dif, ...
    max_deg_pol_x, max_deg_pol_y, max_deg_pol_dif, trans, ...
    correction_y, correction_x, correction_dif)
    
    %% Initialization
    PSF_nb_x = 2*rad_ROI+1 ;
    PSF_nb_y = 2*rad_ROI+1 ;

    % Number of wavelengths
    nb_lambda = length(list_lambda) ;
    
    % Number of lenslets
    nb_lenslet = size(list_x_0, 1) ;
        
    % Transmission
    if nargin < 16 || isempty(trans)
        trans = ones(nb_lenslet, 1) ;
    end
    if size(trans, 2) == 1
        % Duplicating transmission at all wavelength
        trans = repmat(trans, ...
            [1,nb_lambda]) ;
    end
    
    % correction_y
    if nargin < 17 || isempty(correction_y)
        correction_y = zeros(nb_lenslet, deg_pol_y(3)+1) ;
    elseif size(correction_y, 2)<deg_pol_y(3)+1
        correction_y = cat(2, ...
            correction_y, ...
            zeros(nb_lenslet, deg_pol_y(3)+1-size(correction_y, 2))) ;
    end
    
    % correction_x
    if nargin < 18 || isempty(correction_x)
        correction_x = zeros(nb_lenslet, deg_pol_x(3)+1) ;
    elseif size(correction_x, 2)<deg_pol_x(3)+1
        correction_x = cat(2, ...
            correction_x, ...
            zeros(nb_lenslet, deg_pol_x(3)+1-size(correction_x, 2))) ;
    end
    
    % correction_dif
    if nargin < 19 || isempty(correction_dif)
        correction_dif = zeros(nb_lenslet, deg_pol_dif(3)+1) ;
    end

    % Positions of the pixels on the sensor
    sensor_x = get_Fourier_vector(pix.nb_x, pix.dx) ;
    sensor_y = get_Fourier_vector(pix.nb_y, pix.dy) ;

    %% Building of the sparse projection matrix SM
    ind_list_yx = cell(nb_lambda*nb_lenslet, 1) ;
    ind_list_lensletl = cell(nb_lambda*nb_lenslet, 1) ;
    ind_val = cell(nb_lambda*nb_lenslet, 1) ;
    delta_per = 0.5 ;
    [per_aux, lastprint] = display_percentage('init', ...
        '   Building matrix') ;
    % Loop on the lenslet
    for lenslet = 1:nb_lenslet
        %% Percentage
        [per_aux, lastprint] = display_percentage('iter', ...
            {(lenslet-1)/nb_lenslet*100, per_aux, delta_per, ...
            lastprint}) ;

        %% Determination of the positions and diffraction law of the
        % current spectrum
        % y-law
        Pol = LinOpPolynomial( ...
            {list_y_0(lenslet)*ones(nb_lambda,1), ...
            list_x_0(lenslet)*ones(nb_lambda,1), ...
            reshape(list_lambda, [nb_lambda,1])}, ...
            deg_pol_y, max_deg_pol_y) ;
        Pol_y = Pol*coef_pol_y ;
        
        % x-law
        Pol = LinOpPolynomial( ...
            {list_y_0(lenslet)*ones(nb_lambda,1), ...
            list_x_0(lenslet)*ones(nb_lambda,1), ...
            reshape(list_lambda, [nb_lambda,1])}, ...
            deg_pol_x, max_deg_pol_x) ;
        Pol_x = Pol*coef_pol_x ;
        
        % dif-law
        Pol = LinOpPolynomial( ...
            {list_y_0(lenslet)*ones(nb_lambda,1), ...
            list_x_0(lenslet)*ones(nb_lambda,1), ...
            reshape(list_lambda, [nb_lambda,1])}, ...
            deg_pol_dif, max_deg_pol_dif) ;
        Pol_dif = Pol*coef_pol_dif ;
        
        % Local correction
        Pol_y = Pol_y + LinOpPolynomial( ...
            {reshape(list_lambda, [nb_lambda,1])}, ...
            deg_pol_y(3), deg_pol_y(3))*correction_y(lenslet,:)' ;
        Pol_x = Pol_x + LinOpPolynomial( ...
            {reshape(list_lambda, [nb_lambda,1])}, ...
            deg_pol_x(3), deg_pol_x(3))*correction_x(lenslet,:)' ;
        Pol_dif = Pol_dif + LinOpPolynomial( ...
            {reshape(list_lambda, [nb_lambda,1])}, ...
            deg_pol_dif(3), deg_pol_dif(3))*correction_dif(lenslet,:)' ;
        
        %% Loop on the wavelengths
        for l = 1:nb_lambda
            % Position of the current wavelength in the lenslet
            p_x = Pol_x(l) ;
            p_y = Pol_y(l) ;

            % Closest position on the grid
            [~, ind_x] = min(abs(sensor_x-p_x)) ;
            [~, ind_y] = min(abs(sensor_y-p_y)) ;

            % List of indexes to compute the elementary PSF
            ind_x = ind_x+get_Fourier_vector(PSF_nb_x, 1) ;
            ind_y = ind_y+get_Fourier_vector(PSF_nb_y, 1) ;

            % Keeping index on the grid
            ind_x = unique(max(min(ind_x,pix.nb_x),1)) ;
            ind_y = unique(max(min(ind_y,pix.nb_y),1)) ;
            list_x = sensor_x(ind_x)' ;
            list_y = sensor_y(ind_y)' ;
            [ind_xx, ind_yy] = meshgrid(ind_x, ind_y) ;

            % Computation of the elementary PSF
            val = get_pattern_simulation(list_y, list_x, p_y, ...
                p_x, Pol_dif(l), 1, 0, pattern_model) ;

            % Update of the global sparse matrix
            ind_list_yx{lenslet+(l-1)*nb_lenslet} = ...
                ind_yy(:) + pix.nb_y*(ind_xx(:)-1) ;
            ind_list_lensletl{lenslet+(l-1)*nb_lenslet} = ...
                ones(size(ind_yy(:)))* ...
                (lenslet + nb_lenslet*(l-1)) ;
            ind_val{lenslet+(l-1)*nb_lenslet} = ...
                trans(lenslet, l)*val(:) ;
        end
    end

    %% Building finish
    fprintf(repmat('\b', 1, lastprint));
    lastprint = fprintf(['   Building matrix: declaration of', ...
        ' the sparse matrix...']);

    %% Declaration of the sparse projection matrix SM
    ind_list_yx= cell2mat(ind_list_yx) ;
    ind_list_lensletl= cell2mat(ind_list_lensletl) ;
    ind_val= cell2mat(ind_val) ;
    Spec_proj = sparse(ind_list_yx, ind_list_lensletl, ind_val, ...
        pix.nb_y*pix.nb_x, nb_lenslet*nb_lambda) ;
    Spec_proj = LinOpSparseMatrix(Spec_proj, [nb_lenslet, nb_lambda], ...
        [pix.nb_y, pix.nb_x]) ;

    %% Cleaning percentage
    display_percentage('exit', lastprint) ;
    fprintf(repmat('\b', 1, ...
        length('   Building matrix: done!  '))) ;
end