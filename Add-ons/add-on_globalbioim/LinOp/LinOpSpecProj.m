classdef LinOpSpecProj < LinOp
    % LinOpSpecProj: Projection matrix to project the spectra diffracted by
    % the lenslet array on the sensor
    % 
    % :param lenslet_grid: lenslet grid
    %   - lenslet_grid.[list_x, list_y] -> List of the center of the
    %   lenslet ([x,y]-axis) on the sensor in the different calibration
    %   wavelengths
    %   - lenslet_grid.trans -> List of the transmission of each lenslet
    %   and possibly at all wavelength (default = 1)
    %
    % :param sensor: cartesian grid on the sensor
    %   - sensor.[nb_x, nb_y] -> The size of the domain
    %   - sensor.[dx, dy] -> The pixel side pitch
    %
    % :param spectra: spectra features
    %   - spectra.list_lambda -> List of the wavelengths in the model
    %   - spectra.list_lambda_cal -> List of the calibration wavelengths
    %   - spectra.deg_pol_dif -> Degree of the polynomial to fit the
    %   diffraction pattern extension
    %   - spectra.deg_pol_spec=[deg_pol_spec_xl, deg_pol_spec_y] -> Degree
    %   of the polynomial to fit the spectra position along the dispersion
    %   direction(xl) and its orthogonal direction (y) obtained with the
    %   calibration data
    %   - spectra.sig_cal=[sig_xl, sig_y] -> List of the spreading of the
    %   psf along the dispersion direction(xl) and its orthogonal
    %   direction (y) obtained with the calibration data
    %   - spectra.[PSF_nb_x, PSF_nb_y] -> Spatial extension of the 
    %   elementary PSF for a given lenslet and for a given wavelength
    %
    % All attributes of parent class :class:`LinOp` are inherited. 
    %
    % **Example** Spec_proj = get_spec_proj(lenslet_grid, sensor,
    %   lambda)
    %
    % See also :class:`Map` :class:`LinOp`
  
    %%    Copyright (C) 2018
    %     Created: 04/10/2018 (mm/dd/yyyy)
    %     Modified: 04/27/2018 (mm/dd/yyyy) Taking in account that the
    %     dispersion is a list depending of the wavelength.
    %     Modified: 04/27/2018 (mm/dd/yyyy) New definition of the
    %     polynomial fitting
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
        function this = LinOpSpecProj(lenslet_grid, sensor, spectra)
            %% Initialization
            % Number of wavelengths
            nb_lambda = length(spectra.list_lambda) ;
            % Number of wavelengths used in the calibration
            nb_lambda_cal = length(spectra.list_lambda_cal) ;
            % Number of lenslets
            nb_lenslet = size(lenslet_grid.list_x, 1) ;
            % Global orientation of the spectra to take into
            % account the spatial extensions given by sig_cal (deg)
            theta_spec = lenslet_grid.list_x(:,2:nb_lambda_cal) - ...
                lenslet_grid.list_x(:,1) + ...
                1i*(lenslet_grid.list_y(:,2:nb_lambda_cal) - ...
                lenslet_grid.list_y(:,1)) ;
            theta_spec = 180/pi*angle(mean(theta_spec(:)))  ;
            if ~isfield(lenslet_grid, 'trans')
                lenslet_grid.trans = ones(nb_lenslet, 1) ;
            end
            if size(lenslet_grid.trans, 2) == 1
                % Duplicating transmission at all wavelength
                lenslet_grid.trans = repmat(lenslet_grid.trans, ...
                    [1,nb_lambda]) ;
            end
            
            % Parameters of this
            this.name = 'LinOpSpecProj' ;
            this.isInvertible = false ;
            this.sizein = [nb_lenslet, nb_lambda] ;
            this.sizeout = [sensor.nb_y, sensor.nb_x] ;
            
            % Simulation parameters for the diffraction patterns
            sig_fit = fit_pol(spectra.deg_pol_dif, ...
                spectra.list_lambda_cal, spectra.sig_cal, ...
                spectra.list_lambda) ;
            
            % Positions of the pixels on the sensor / Extension of the
            % domain to properly normalize the elementary PSFs which fall
            % on the edges
            nb_x_extended = sensor.nb_x+2*spectra.PSF_nb_x ;
            nb_y_extended = sensor.nb_y+2*spectra.PSF_nb_y ;
            sensor_x = get_Fourier_vector(nb_x_extended, sensor.dx) ;
            sensor_y = get_Fourier_vector(nb_y_extended, sensor.dy) ;
                        
            %% Building of the sparse projection matrix SM
            ind_list_yx = cell(nb_lambda*nb_lenslet, 1) ;
            ind_list_lensletl = cell(nb_lambda*nb_lenslet, 1) ;
            ind_val = cell(nb_lambda*nb_lenslet, 1) ;
            lastprint = fprintf('   Building matrix...') ;
            delta_per = 0.5 ;
            per_aux = 0 ;
            % Loop on the lenslet
            for lenslet = 1:nb_lenslet
                %% Percentage
                per = (lenslet-1)/nb_lenslet*100 ;
                if per>per_aux
                    fprintf(repmat('\b', 1, lastprint));
                    lastprint = fprintf(['   Building matrix: ', ...
                        num2str(per_aux), ' %%']);
                    per_aux = (floor(per/delta_per)+1)*delta_per ;
                end
                
                %% Determination of the positions of the current spectrum
                % Calibration spectrum
                list_x = lenslet_grid.list_x(lenslet,:)' ;
                list_y = lenslet_grid.list_y(lenslet,:)' ;
                
                % Fitted spectrum
                [list_x, list_y] = fit_spectrum(spectra.deg_pol_spec, ...
                    [list_x, list_y], spectra.list_lambda_cal, ...
                    spectra.list_lambda) ;
                
                %% Loop on the wavelengths
                for l = 1:nb_lambda
                    % Position of the current wavelength in the lenslet
                    p_x = list_x(l) ;
                    p_y = list_y(l) ;
                    
                    % Closest position on the grid
                    [~, ind_x] = min(abs(sensor_x-p_x)) ;
                    [~, ind_y] = min(abs(sensor_y-p_y)) ;
                    
                    % List of indexes to compute the elementary PSF
                    ind_x = ind_x+get_Fourier_vector(spectra.PSF_nb_x, 1) ;
                    ind_y = ind_y+get_Fourier_vector(spectra.PSF_nb_y, 1) ;
                    
                    % Keeping index on the grid
                    ind_x = unique(max(min(ind_x,nb_x_extended),1)) ;
                    ind_y = unique(max(min(ind_y,nb_y_extended),1)) ;
                    
                    % Rotation in the frame of the slit
                    [ind_xx, ind_yy] = meshgrid(ind_x, ind_y) ;
                    xx = sensor_x(ind_xx)-p_x ;
                    yy = sensor_y(ind_yy)-p_y ;
                    [xx, yy] = rot_2D(-theta_spec, xx(:), yy(:)) ;
                    
                    % Computation of the elementary PSF
                    val = exp(-0.5*((xx/sig_fit(l,1)).^2 + ...
                        (yy/sig_fit(l,2)).^2)) ;
                    
                    % Normalization
                    if sum(val(:))>0
                        val = val./sum(val(:)) ;
                    end
                    
                    % Cropping back to the initial domain
                    ind_yy = ind_yy(:)-spectra.PSF_nb_y ;
                    sel_ind = (ind_yy>0) & (ind_yy<sensor.nb_y+1) ;
                    ind_xx = ind_xx(:)-spectra.PSF_nb_x;
                    sel_ind = sel_ind & ...
                        (ind_xx>0) & (ind_xx<sensor.nb_x+1) ;
                    
                    % Update of the global sparse matrix
                    ind_list_yx{lenslet+(l-1)*nb_lenslet} = ...
                        ind_yy(sel_ind) + sensor.nb_y*(ind_xx(sel_ind)-1) ;
                    ind_list_lensletl{lenslet+(l-1)*nb_lenslet} = ...
                        ones(size(ind_yy(sel_ind)))* ...
                        (lenslet + nb_lenslet*(l-1)) ;
                    ind_val{lenslet+(l-1)*nb_lenslet} = ...
                        lenslet_grid.trans(lenslet, l)*val(sel_ind) ;
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
            this.SM = sparse(ind_list_yx, ind_list_lensletl, ind_val, ...
                sensor.nb_y*sensor.nb_x, nb_lenslet*nb_lambda) ;
            
            %% Cleaning percentage
            fprintf(repmat('\b', 1, lastprint)) ;
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