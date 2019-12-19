classdef OpL2norm < Map
    % OpL2norm: Operator computing the squared L2-norm of a vector.
    % 
    % :param sizein: size of the input matrix
    % 
    % :param index: index along which the norm is computed (default: the
    %   last one)
    %
    % **Example** L2norm = OpL2norm(sizein, index)
    %
    % **Example** L2norm = OpL2norm(sizein)
    %
    % See also :class:`Map`
  
    %%    Copyright (C) 2018
    %     Created: 05/03/2018 (mm/dd/yyyy)
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
        index ; % index along which the norm is computed
    end
    
    %% Constructor
    methods
        function this = OpL2norm(sizein, index)
            this.name = 'OpL2norm' ;
            this.isDifferentiable=true;
            this.isInvertible = false ;
            
            
            this.index = length(sizein) ;
            if nargin > 1
                if ~isempty(index)
                    this.index = index ;
                end
            end
            this.sizein = sizein ;
            this.sizeout = sizein ;
            this.sizeout(this.index) = [] ;
            if length(this.sizeout) == 1 % The output is a vector
                if this.index==1 % The output is a row vector
                    this.sizeout = [1, this.sizeout] ;
                else % The output is a column vector
                    this.sizeout = [this.sizeout, 1] ;
                end
            end
        end
    end
	
    %% Core Methods containing implementations (Protected)
	methods (Access = protected)
        function y = apply_(this,x)
            % Reimplemented from parent class :class:`Map`.
            
            y = sum(x.^2,this.index) ;
        end
        
        function x = applyJacobianT_(this,y,v)
            % Reimplemented from parent class :class:`Map`.
            
            % Permutation of the projection dimension to the first
            % dimention
            ind_perm = [this.index, setdiff(1:length(this.sizein), ...
                this.index)] ;
            v = permute(v, ind_perm) ;
            
            % Loop on the projection dimension
            if this.sizein(this.index)>1
                sizeaux = [this.sizein(this.index), this.sizeout] ;
                if this.sizeout(1) == 1 % Row vector
                    sizeaux(2) = [] ;
                end
                x = zeros(sizeaux) ;
            else
                sizeaux = this.sizein ;
                sizeaux(this.index) = [] ;
                x = zeros([1, sizeaux]) ;
            end
            for i = 1:this.sizein(this.index)
                x(i,:) = 2*v(i,:).*y(:)' ;
            end
            % Inverse permutation to go back to the initial order
            x = ipermute(x, ind_perm) ;
        end
    end
end
