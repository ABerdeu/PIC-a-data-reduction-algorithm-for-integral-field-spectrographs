classdef LinOpSparseMatrix < LinOp
    % LinOpSparseMatrix: Equivalent to the LinOpMatrix for sparse matrix
    % 
    % :param SM: the sparse matrix (which is 2D in Matlab)
    %
    % :param sizein: the size of the input of the sparse matrix
    %
    % :param sizeout: the size of the output of the sparse matrix
    %
    % All attributes of parent class :class:`LinOp` are inherited. 
    %
    % **Example** Sparse_Mat = LinOpSparseMatrix(SM, sizein, sizeout)
    %
    % See also :class:`Map` :class:`LinOp`
  
    %%    Copyright (C) 2018
    %     Created: 11/19/2018 (mm/dd/yyyy)
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
        SM ;        % Sparse matrix
    end
    
    %% Constructor
    methods
        function this = LinOpSparseMatrix(SM, sizein, sizeout)
            
            % Parameters of this
            this.name = 'LinOpSparseMatrix' ;
            this.isInvertible = false ;
            this.sizein = sizein ;
            this.sizeout = sizeout ;
            this.SM = SM ;
            
            % Checking sizes
            if size(SM, 1)~=prod(sizeout)
                error(['The total size of the output (', ...
                    num2str(prod(sizeout)), ...
                    ') and the number of rows of the matrix (', ...
                    num2str(size(SM, 1)), ...
                    ') do not match...']) ;
            end
            if size(SM, 2)~=prod(sizein)
                error(['The total size of the input (', ...
                    num2str(prod(sizein)), ...
                    ') and the number of columns of the matrix (', ...
                    num2str(size(SM, 2)), ...
                    ') do not match...']) ;
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