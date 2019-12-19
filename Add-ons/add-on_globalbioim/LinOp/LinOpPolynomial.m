classdef LinOpPolynomial < LinOp
    % LinOpPolynomial: Linear operator on the polynomial space
    % 
    % :param list_coord: cell array containing the list of the
    % coordinates on which the polynomial is estimated. If there is only
    % coordinate, it is not necessary a cell array.
    %
    % :param degree: degree of the polynomial expression on each dimension.
    % The maximal degree of the polynomial expression is then sum(degree)
    % If only one value is given, the polynomial expression is assumed to
    % be of total degree 'degree'.
    %
    % :param max_degree: maximal total degree of the polynomial (default:
    % the value is given by the list of 'degree') ;
    %
    % All attributes of parent class :class:`LinOp` are inherited. 
    %
    % **Note**: it is implemented for dimensions lower than 3
    %
    % **Example** Interp = LinOpPolynomial(sizeout, degree)
    %
    % See also :class:`Map` :class:`LinOp`
  
    %%    Copyright (C) 2018
    %     Created: 10/08/2018 (mm/dd/yyyy)
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
        list_elemPol ;  % list of the estimation on the grid of the
            % elementary polynomials [1, X, X^2, ...]
        list_deg ;      % list of the degree of each coefficient
    end
    
    %% Constructor
    methods
        function this = LinOpPolynomial(list_coord, degree, max_degree)
            %% Initialization
            this.name = 'LinOpPolynomial' ;             
            this.isInvertible = false ;
            
            % Checking input
            if ~isa(list_coord, 'cell')
                list_coord = {list_coord} ;
            end
            this.sizeout = size(list_coord{1}) ;
            ndms = length(list_coord) ;
            
            % degree of the polynomial
            if nargin < 3 || isempty(max_degree)
                if length(degree) == 1
                    max_degree = degree ;
                    degree = repmat(degree, ndms) ;
                else
                    max_degree = sum(degree) ;
                end
            end
            
            % List of degrees for the elementary polynomials
            switch ndms
                case 1 % 1D grid
                    list_deg = (0:degree)' ;
                
                case 2 % 2D grid
                    [list_deg_2, list_deg_1] = meshgrid( ...
                        0:degree(2), ...
                        0:degree(1)) ;
                    list_deg = [list_deg_1(:), list_deg_2(:)] ;
                    list_deg = list_deg(sum(list_deg,2)<(max_degree+1),:) ;
                
                case 3 % 3D grid
                    [list_deg_2, list_deg_1, list_deg_3] = meshgrid( ...
                        0:degree(2), ...
                        0:degree(1), ...
                        0:degree(3)) ;
                    list_deg = [list_deg_1(:), list_deg_2(:), ...
                        list_deg_3(:)] ;
                    list_deg = list_deg(sum(list_deg,2)<(max_degree+1),:) ;
                    
                otherwise
                    error(['The number of dimensions is ', ...
                        num2str(ndms), ...
                        '... It must be lower than 3...']) ;
            end
            
            % Number of coefficient of the polynomial
            this.list_deg = list_deg ;
            this.sizein = [size(list_deg, 1), 1] ;
            
            % Defining the elementary polynomials
            this.list_elemPol = cell(this.sizein(1), 1) ;
            switch ndms
                case 1 % 1D grid
                    % Coordinates on the grid
                    x = list_coord{1} ;
                    
                    % Elementary polynomials
                    for d = 1:this.sizein
                        this.list_elemPol{d} = x.^list_deg(d) ;
                    end
                    
                case 2 % 2D grid
                    % Coordinates on the grid
                    y = list_coord{1} ;
                    x = list_coord{2} ;
                    if ~isequal(size(y), size(x))
                        error(['The coordinates must have the same', ...
                            ' size...']) ;
                    end
                    
                    % Elementary polynomials
                    for d = 1:this.sizein
                        this.list_elemPol{d} = x.^list_deg(d,1) .* ...
                            y.^list_deg(d,2) ;
                    end
                
                case 3 % 3D grid
                    % Coordinates on the grid
                    y = list_coord{1} ;
                    x = list_coord{2} ;
                    z = list_coord{3} ;
                    if ~isequal(size(y), size(x)) || ...
                            ~isequal(size(y), size(z))
                        error(['The coordinates must have the same', ...
                            ' size...']) ;
                    end
                    
                    % Elementary polynomials
                    for d = 1:this.sizein
                        this.list_elemPol{d} = x.^list_deg(d,1) .* ...
                            y.^list_deg(d,2) .* ...
                            z.^list_deg(d,3) ;
                    end
            end
        end
    end
	
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyAdjoint_(this,x)
	methods (Access = protected)
        %% apply_
        function y = apply_(this,x)
            % Initialization
            y = zeros(this.sizeout) ;
            
            % Loop on the elementary polynomials
            for d = 1:this.sizein
                y = y + x(d)*this.list_elemPol{d} ;
            end
        end
        
        %% applyAdjoint_
        function y = applyAdjoint_(this,x)
            % Initialization
            y = zeros(this.sizein) ;
            
            % Loop on the elementary polynomials
            for d = 1:this.sizein
                y_d = x.*this.list_elemPol{d} ;
                y(d) = sum(y_d(:)) ;
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