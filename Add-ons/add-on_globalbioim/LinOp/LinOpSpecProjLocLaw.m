classdef LinOpSpecProjLocLaw < LinOp
    % LinOpSpecProjLocLaw: Projection matrix to project the spectra
    % diffracted by the lenslet array on the sensor based on local laws for
    % the spectra projection (gaussian patterns)
    % 
    % :param sensor: cartesian grid on the sensor
    %   - sensor.[nb_x, nb_y] -> The size of the domain
    %   - sensor.[dx, dy] -> The pixel side pitch
    %
    % :param spectra: spectra features
    %   - spectra.list_lambda -> List of the wavelengths in the model
    %   - spectra.[PSF_nb_x, PSF_nb_y] -> Spatial extension of the 
    %   elementary PSF for a given lenslet and for a given wavelength
    %   - spectra.[coef_pol_y, coef_pol_x] -> Coefficients of the
    %   polynomial laws ([x,y]-axis) on the sensor
    %   - spectra.coef_pol_dif -> Coefficients of the diffraction law
    %   of the different spectra
    %   - spectra.trans -> List of the transmission of each lenslet
    %   and possibly at all wavelength (default = 1)
    %
    % :pattern_model: flags to specify the pattern to fit
    %   pattern_model.flag_norm -> Normalized?
    %   pattern_model.flag_profile -> Gaussian / Hexagon
    %   pattern_model.oversampling -> oversampling of the pattern 
    %   (default: 1)
    %
    % :flag_verbose: flag to display the percentage (default: true)
    %
    % All attributes of parent class :class:`LinOp` are inherited. 
    %
    % **Example** Spec_proj = get_spec_proj(sensor, spectra)
    %
    % See also :class:`Map` :class:`LinOp`
  
    %%    Copyright (C) 2018
    %     Created: 04/10/2018 (mm/dd/yyyy)
    %     Modified: 04/27/2018 (mm/dd/yyyy) Taking in account that the
    %     dispersion is a list depending of the wavelength.
    %     Modified: 04/27/2018 (mm/dd/yyyy) New definition of the
    %     polynomial fitting
    %     Modified: 09/06/2018 (mm/dd/yyyy) The projection pattern is
    %     axisymmetric
    %     Modified: 03/11/2018 (mm/dd/yyyy) Local law
    %     Anthony Berdeu (Laboratoire Hubert Curien)
    %
    %     This program is free software: you can redistribute it and/or modify
    %     it under the terms of the GNU General Public License as published by
    %     the Free Software Foundation, either version 3 of the License, or
    %     (at your option) any later version.
    %
    %     This program is distributed in the hope that it will be useful,
    %     but WITHOUT ANY WARRANTY; without even the implied warranty of
    %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %     GNU General Public License for more details.
    %
    %     You should have received a copy of the GNU General Public License
    %     along with this program.  If not, see <http://www.gnu.org/licenses/>.

    properties
        SM ;        % Sparse matrix of the spectra projections
    end
    
    %% Constructor
    methods
        function this = LinOpSpecProjLocLaw(sensor, spectra, ...
            pattern_model, flag_verbose)
            %% Initialization
            if nargin < 4
                flag_verbose = true ;
            end
            % List of wavelength
            list_lambda = spectra.list_lambda ;
            
            % Number of wavelengths
            nb_lambda = length(list_lambda) ;
            
            % Coefficients of the different polynomial laws
            coef_pol_y = spectra.coef_pol_y ;
            coef_pol_x = spectra.coef_pol_x ;
            coef_pol_dif = spectra.coef_pol_dif ;
            
            % Number of coefficients
            deg_pol_y = size(coef_pol_y, 2) - 1 ;
            deg_pol_x = size(coef_pol_x, 2) - 1 ;
            deg_pol_dif = size(coef_pol_dif, 2) - 1 ;
            
            % y-law polynomial
            Pol_y = zeros(nb_lambda, deg_pol_y+1) ;
            for d = 0:deg_pol_y
                Pol_y(:,d+1) = list_lambda.^d ;
            end

            % x-law polynomial
            Pol_x = zeros(nb_lambda, deg_pol_x+1) ;
            for d = 0:deg_pol_x
                Pol_x(:,d+1) = list_lambda.^d ;
            end

            % Diffraction polynomial
            Pol_dif = zeros(nb_lambda, deg_pol_dif+1) ;
            for d = 0:deg_pol_dif
                Pol_dif(:,d+1) = list_lambda.^d ;
            end
            
            % Number of lenslets
            nb_lenslet = size(coef_pol_x, 1) ;
            
            if ~isfield(spectra, 'trans')
                trans = ones(nb_lenslet, 1) ;
            else
                trans = spectra.trans ;
            end
            if size(trans, 2) == 1
                % Duplicating transmission at all wavelength
                trans = repmat(trans, ...
                    [1,nb_lambda]) ;
            end
            
            % Parameters of this
            this.name = 'LinOpSpecProjLocLaw' ;
            this.isInvertible = false ;
            this.sizein = [nb_lenslet, nb_lambda] ;
            this.sizeout = [sensor.nb_y, sensor.nb_x] ;
            
            % Positions of the pixels on the sensor
            sensor_x = get_Fourier_vector(sensor.nb_x, sensor.dx) ;
            sensor_y = get_Fourier_vector(sensor.nb_y, sensor.dy) ;
                        
            %% Building of the sparse projection matrix SM
            ind_list_yx = cell(nb_lambda*nb_lenslet, 1) ;
            ind_list_lensletl = cell(nb_lambda*nb_lenslet, 1) ;
            ind_val = cell(nb_lambda*nb_lenslet, 1) ;
            delta_per = 0.5 ;
            if flag_verbose
                [per_aux, lastprint] = display_percentage('init', ...
                    '   Building matrix') ;
            end
            % Loop on the lenslet
            for lenslet = 1:nb_lenslet
                %% Percentage
                if flag_verbose
                    [per_aux, lastprint] = display_percentage('iter', ...
                        {(lenslet-1)/nb_lenslet*100, per_aux, ...
                        delta_per, lastprint}) ;
                end
                
                %% Determination of the positions of the current spectrum
                % Calibration spectrum
                list_y_spec = Pol_y*coef_pol_y(lenslet,:)' ;
                list_x_spec = Pol_x*coef_pol_x(lenslet,:)' ;
                list_dif_spec = Pol_dif*coef_pol_dif(lenslet,:)' ;
                
                %% Loop on the wavelengths
                for l = 1:nb_lambda
                    % Position of the current wavelength in the lenslet
                    p_x = list_x_spec(l) ;
                    p_y = list_y_spec(l) ;
                    
                    % Closest position on the grid
                    [~, ind_x] = min(abs(sensor_x-p_x)) ;
                    [~, ind_y] = min(abs(sensor_y-p_y)) ;
                    
                    % List of indexes to compute the elementary PSF
                    ind_x = ind_x+get_Fourier_vector(spectra.PSF_nb_x, 1) ;
                    ind_y = ind_y+get_Fourier_vector(spectra.PSF_nb_y, 1) ;
                    
                    % Keeping index on the grid
                    ind_x = unique(max(min(ind_x,sensor.nb_x),1)) ;
                    ind_y = unique(max(min(ind_y,sensor.nb_y),1)) ;
                    list_x = sensor_x(ind_x)' ;
                    list_y = sensor_y(ind_y)' ;
                    [ind_xx, ind_yy] = meshgrid(ind_x, ind_y) ;
                    
                    % Computation of the elementary PSF
                    val = get_pattern_simulation(list_y, list_x, p_y, ...
                        p_x, list_dif_spec(l), 1, 0, pattern_model) ;
                    
                    % Update of the global sparse matrix
                    ind_list_yx{lenslet+(l-1)*nb_lenslet} = ...
                        ind_yy(:) + sensor.nb_y*(ind_xx(:)-1) ;
                    ind_list_lensletl{lenslet+(l-1)*nb_lenslet} = ...
                        ones(size(ind_yy(:)))* ...
                        (lenslet + nb_lenslet*(l-1)) ;
                    ind_val{lenslet+(l-1)*nb_lenslet} = ...
                        trans(lenslet, l)*val(:) ;
                end
            end
            
            %% Building finish
            if flag_verbose
                fprintf(repmat('\b', 1, lastprint));
                lastprint = fprintf(['   Building matrix: ', ...
                    'declaration of the sparse matrix...']);
            end
            
            %% Declaration of the sparse projection matrix SM
            ind_list_yx= cell2mat(ind_list_yx) ;
            ind_list_lensletl= cell2mat(ind_list_lensletl) ;
            ind_val= cell2mat(ind_val) ;
            this.SM = sparse(ind_list_yx, ind_list_lensletl, ind_val, ...
                sensor.nb_y*sensor.nb_x, nb_lenslet*nb_lambda) ;
            
            %% Cleaning percentage
            if flag_verbose
                display_percentage('exit', lastprint) ;
                fprintf(repmat('\b', 1, ...
                    length('   Building matrix: done!  '))) ;
            end
        end
    end
	
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyAdjoint_(this,x)
	methods (Access = protected)
        %% apply_
        function y = apply_(this,x)   
            % Reimplemented from parent class :class:`LinOp`.
            y = this.SM*x(:) ;
            y = reshape(y, this.sizeout) ;
        end
        
        %% applyAdjoint_
        function y = applyAdjoint_(this,x)            
            % Reimplemented from parent class :class:`LinOp`.
            y = this.SM'*x(:) ;
            y = reshape(y, this.sizein) ;
        end
        
        %% applyHtH_
%         function y = applyHtH_(this,x)
%             % Reimplemented from parent class :class:`LinOp`.
%         end
        
        %% applyHHt_
%         function y = applyHHt_(this,x)            
%             % Reimplemented from parent class :class:`LinOp`.
%         end
		
        %% makeHtH_
%         function M = makeHtH_(this)
%             % Reimplemented from parent class :class:`LinOp`.
%         end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get the coordinates in the Fourier space.
% Because of the fft and the fftshift properties, a zero must be in the
% list at a specific position.
% 
% Created: 03/19/2015 (mm/dd/yyyy)
% Modified: 10/26/2015 (mm/dd/yyyy) (changing the name, from
% "get_kernel_vector")
% Author: Anthony Berdeu (PhD student at CEA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #nb_el# the number of elements in the vector
%
% #d_el# the elementary length of the elements in the vector
%
%%%%%%%%
% Ouput
%%%%%%%%
% #pos_Fourier# the position of the kernel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pos_Fourier] = get_Fourier_vector(nb_el, d_el)

    if mod(nb_el,2)
        % There is an odd number of elements, the center of the frame is at
        % the center of the vector
        pos_Fourier = ((0:nb_el-1)-(nb_el-1)/2)*d_el ;
    else
        % There is an even number of elements, the center of the frame is
        % at the position round(nb_pix/2)+1
        pos_Fourier = ((0:nb_el-1)-nb_el/2)*d_el ;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to rotate a pair of matrix.
% 
% Created: 04/10/2018 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #theta# the rotation angle (deg)
%
% #[x_in, y_in]# the input vectors to rotate of the angle theta
%
%%%%%%%%
% Ouput
%%%%%%%%
% #[x_out, y_out]# the output vectors rotated of the angle theta
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_out, y_out] = rot_2D(theta, x_in, y_in)
    x_out = cosd(theta).*x_in - sind(theta).*y_in ;
    y_out = sind(theta).*x_in + cosd(theta).*y_in ;
end