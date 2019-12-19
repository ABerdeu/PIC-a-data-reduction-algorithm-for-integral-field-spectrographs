classdef LinOpInterp < LinOp
    % LinOpInterp: Linear interpolation operator on a cartesian grid
    % 
    % :param coord_in: coordinates of the sample points.
    %
    % :param coord_out: coordinates of the query points.
    %
    % :param extrapVal: extrapolation value for the query points which fall
    % outside the intial grid (default = 0, see notes).
    %
    % :param sizein: dimension of the hyper volume (default = dimensions
    % given by coord_in)
    %
    % :param interp_dim: dimensions of the hyper surface along which the
    % interpolation must be performed (default ! the first dimensions of
    % the hyper volume)
    %
    % All attributes of parent class :class:`LinOp` are inherited. 
    %
    % **Note**: it is possible to extrapolate with values different from 0
    % but in that case the interpolation operation is not linear.
    %
    % **Example** Interp = LinOpInterp(coord_in, coord_out)
    %
    % **Example** Interp = LinOpInterp(coord_in, coord_out, extrapVal)
    %
    % **Example** Interp = LinOpInterp(coord_in, coord_out, extrapVal, sizein)
    %
    % **Example** Interp = LinOpInterp(coord_in, coord_out, extrapVal, sizein, interp_dim)
    %
    % **Example** Interp = LinOpInterp(coord_in, coord_out, [], sizein)
    %
    % **Example** Interp = LinOpInterp(coord_in, coord_out, [], sizein, interp_dim)
    %
    % See also :class:`Map` :class:`LinOp`
  
    %%    Copyright (C) 2018
    %     Created: 03/27/2018 (mm/dd/yyyy)
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
        ndms ;      % number of dimensions of the hyper surface
        ndms_vol ;  % number of dimensions of the hyper volume
        interp_dim ;% interpolation dimensions
        weight ;    % weighting coefficients of the linear interpolation
        V_ind       % index of the interpolated values
        inner_ind   % index of the query positions which lie inside the
            % grid
        extrapVal ; % extrapolation value for the query points which fall
            % outside the intial grid.
    end
    
    %% Constructor
    methods
        function this = LinOpInterp(varargin)
            %% Initialization
            this.name = 'LinOpInterp' ;             
            this.isInvertible = false ;
            
            % Coordinates of the sample points.
            coord_in = varargin{1} ;
            this.ndms = length(coord_in) ;
            
            % Determination of sizein
            if nargin>3
                this.sizein = varargin{4} ;
            else
                this.sizein = zeros(1, this.ndms) ;
                for dim = 1:this.ndms
                    this.sizein(dim) = length(coord_in{dim}) ;
                end
            end
            this.ndms_vol = length(this.sizein) ;
            
            % Interpolation dimensions
            if nargin>4
                this.interp_dim = varargin{5} ;
            else
                this.interp_dim = 1:this.ndms ;
            end
            this.ndms_vol = length(this.sizein) ;
            
            % Query points
            coord_out = varargin{2} ;
            nb_p = size(coord_out,1) ;
            extra_dim = setdiff(1:this.ndms_vol,this.interp_dim) ;
            if isempty(extra_dim)
                this.sizeout = [nb_p, 1] ;
            else
                this.sizeout = [nb_p, this.sizein(extra_dim)] ;
            end
            
            
            
            % Extrapolation value for the query points which fall
            % outside the intial grid.
            if nargin>2 && ~isempty(varargin{3})
                this.extrapVal = varargin{3} ;
            else
                this.extrapVal = 0 ;
            end
            if this.extrapVal ~= 0
                warning(['Extrapolation value different from 0. ', ...
                    'The operator is not linear...']) ;
            end
            
            % Test the number of dimensions
            if size(coord_out,2)~=this.ndms
                error(['Dimensions of the values and the ', ...
                    'coordinates to interpolate do not match!']) ;
            end
            if size(coord_out,2)~=length(this.interp_dim)
                error(['Numbers interpolation dimensions and' ...
                    ' the coordinates to interpolate do not match!']) ;
            end
            
            %% Determination of the interpolated indexes and weighting
            % coefficients
            
            % Initialization
            this.weight = ones(nb_p, this.ndms, 2) ;
            this.V_ind = ones(nb_p, 1, 'uint32') ;
            this.inner_ind = uint32(1:nb_p) ;
            
            % Loop on the query points
            size_aux = [1, this.sizein(this.interp_dim)] ;
            for p = 1:nb_p
                dim = 0 ; % Current dimension
                
                % Loop on the dimensions while the point is considered to
                % be in the grid to build the weighting coefficients and
                % find the index in the value grid
                while dim<this.ndms && this.inner_ind(p)
                    % Increment dim
                    dim = dim+1 ;
                    
                    % Differential positions
                    dp = coord_in{dim}-coord_out(p,dim) ;
                    
                    % index of the first positive value
                    ind_p = find(dp>=0,1) ;

                    % Is the query point lying in the grid on the
                        % p-axis?
                    if dp(1) > 0 ||  isempty(ind_p) % This query point
                            % falls outside the grid
                        this.inner_ind(p) = 0 ;
                        
                    elseif ind_p > 1 % This query point is in the grid
                        ind_m = ind_p-1 ;
                    else % (ind_p==1) && (dp(ind_p)==0) -> point on
                            % the left edge
                        ind_p = 2 ;
                        ind_m = 1 ;
                    end

                    % Determination of the interpolated indexes and
                    % weighting coefficients
                    if this.inner_ind(p)
                        % weight
                        w_p = dp(ind_p)/(dp(ind_p)-dp(ind_m)) ;
                        this.weight(p, dim, :) = [w_p, 1-w_p] ; % p_0 p_1
                            
                        % Update of ind
                        this.V_ind(p) = this.V_ind(p) + ...
                            (ind_m-1)*prod(size_aux(1:dim)) ;
                    end
                end
            end
            
            % Exclusion of the points lying outside the grid
            this.inner_ind = this.inner_ind(this.inner_ind>0) ;
            this.weight = this.weight(this.inner_ind,:,:) ;
            this.V_ind = this.V_ind(this.inner_ind) ;
        end
    end
	
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyAdjoint_(this,x)
	methods (Access = protected)
        %% apply_
        function y = apply_(this,x)
            % Reshape x to put the interpolation dimensions at the
            % beginning
            extra_dim = setdiff(1:this.ndms_vol,this.interp_dim) ;
            if isempty(extra_dim)
                if this.ndms > 1
                    x = permute(x, this.interp_dim) ;
                end
                prod_extra_dim = 1 ;
            else
                x = permute(x, [this.interp_dim, extra_dim]) ;
                prod_extra_dim = prod(this.sizein(extra_dim)) ;
            end
            x = reshape(x, [this.sizein(this.interp_dim), prod_extra_dim]) ;
            size_aux = [1, this.sizein(this.interp_dim)] ;
            nb_ind_surf = prod(size_aux) ;
            
            % Initialization
            y = this.extrapVal*ones(this.sizeout) ;
            if this.extrapVal~=0
                for extra_dim_ind = 1:prod_extra_dim
                    y(this.inner_ind,extra_dim_ind) = 0 ;
                end
            end
            
            % Loop on the dimensions
            ind_aux = zeros(this.ndms, 1) ;
            for dim = 1:this.ndms
                ind_aux(dim) = prod(size_aux(1:dim)) ;
            end
            
            % Loop on the corners of the element
            for corner = 0:(2^this.ndms-1)
                % Initialization
                weight_c = ones(length(this.inner_ind),1) ;
                corners = int2bin(corner, this.ndms) ;
                
                % Loop on the dimensions
                for dim = 1:this.ndms
                    weight_c = weight_c .* ...
                        this.weight(:, dim, corners(dim)+1) ;
                end
                
                % Update of y(this.inner_ind)
                for extra_dim_ind = 1:prod_extra_dim
                    y(this.inner_ind,extra_dim_ind) = ...
                        y(this.inner_ind,extra_dim_ind) + ...
                        weight_c.*x( ...
                        this.V_ind(:) + nb_ind_surf*(extra_dim_ind-1) + ...
                        sum(corners.*ind_aux) ...
                        ) ;
                end
            end
        end
        
        %% applyAdjoint_
        function y = applyAdjoint_(this,x)
            % Initialization
            x_inner = x(this.inner_ind,:) ;
            nb_inner = length(this.inner_ind) ;
            size_aux = [1, this.sizein(this.interp_dim)] ;
            extra_dim = setdiff(1:this.ndms_vol,this.interp_dim) ;
            prod_interp_dim = prod(this.sizein(this.interp_dim)) ;
            if isempty(extra_dim)
                y = zeros([prod_interp_dim, 1]) ;
                prod_extra_dim = 1 ;
            else
                y = zeros([prod_interp_dim, this.sizein(extra_dim)]) ;
                prod_extra_dim = prod(this.sizein(extra_dim)) ;
            end
            
            % Loop on the dimensions
            ind_aux = zeros(this.ndms, 1) ;
            for dim = 1:this.ndms
                ind_aux(dim) = prod(size_aux(1:dim)) ;
            end
            
            for extra_dim_ind = 1:prod_extra_dim
                % List of the indexes and values to accumulate
                list_ind = zeros(nb_inner, 2^this.ndms) ;
                list_V = zeros(nb_inner, 2^this.ndms) ;
                % Loop on the corners of the element
                for corner = 0:(2^this.ndms-1)
                    % Initialization
                    weight_c = ones(length(this.inner_ind),1) ;
                    corners = int2bin(corner, this.ndms) ;

                    % Loop on the dimensions
                    for dim = 1:this.ndms
                        weight_c = weight_c .* ...
                            this.weight(:, dim, corners(dim)+1) ;
                    end

                    % Update of y(this.inner_ind)
                    list_ind(:, corner+1) = this.V_ind(:) + ...
                        sum(corners.*ind_aux) ;
                    list_V(:, corner+1) = weight_c.* ...
                        x_inner(:,extra_dim_ind) ;
                end

                % Adding all contributions
                y(:,extra_dim_ind) = accumarray(list_ind(:), list_V(:), ...
                    [prod_interp_dim, 1]) ;
            end
            
            
            % Reshape y to put the interpolation dimensions at their
            % correct position
            if isempty(extra_dim)
                if this.ndms_vol > 1
                    y = ipermute(y, this.interp_dim) ;
                    y = reshape(y, this.sizein) ;
                end
            else
                y = reshape(y, [this.sizein(this.interp_dim), ...
                    this.sizein(extra_dim)]) ;
                y = ipermute(y, [this.interp_dim, extra_dim]) ;
            end
        end
        
        %% applyHtH_
%         function y = applyHtH_(~,x)
%             % Reimplemented from parent class :class:`LinOp`.
%             y = x ;
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
% Function to convert a positive integer to a binary number
%
% Created: 04/06/2018 (mm/dd/yyyy)
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #m# a positive integer
%
% #n_in# the wanted length of the binary vector (optional)
%
%%%%%%%%
% Ouput
%%%%%%%%
% #b# m written in a binary mode
%
% #n_out# the output length of the binary vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [b, n_out] = int2bin(m, n_in)

    %% Initialization
    n_out = floor(log(m)/log(2))+1 ;
    if nargin>1
        n_out = max(n_out, n_in) ;
    end

    b = zeros(n_out,1) ;
    r = m ;
    
    %% Loop on the powers of 2
    for p = n_out-1:-1:0
        q = floor(r/2^p) ;
        r = r-q*2^p ;
        b(p+1) = q ;
    end
end