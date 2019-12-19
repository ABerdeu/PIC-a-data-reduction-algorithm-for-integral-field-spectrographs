classdef LinOpR2C < LinOp
    % LinOpR2C: To convert a real matrix to a complex matrix. The last
    % dimension of the input matrix is used to store the real and imaginary
    % parts of the final matrix.
    % 
    % :param sz_out: size of the complex output 
    % 
    % All attributes of parent class :class:`LinOp` are inherited. 
    %
    % **Example** R2C = LinOpR2C(sz_out)
    %
    % See also :class:`Map` :class:`LinOp`
  
    %%    Copyright (C) 2018
    %     Created: 05/02/2018 (mm/dd/yyyy)
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

    
    %% Constructor
    methods
        function this = LinOpR2C(sz_out)
            this.name = 'LinOpR2C' ;
            this.isInvertible = true ;
            this.sizeout = sz_out ;
            this.sizein = [sz_out, 2] ;
		end
    end
	
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyAdjoint_(this,x)
	methods (Access = protected)
        function y = apply_(this,x)   
            % Reimplemented from parent class :class:`LinOp`.
            x = reshape(x, [prod(this.sizeout), 2]) ;
            y = x(:,1)+1i*x(:,2) ;
            y = reshape(y, this.sizeout) ;
        end		
        function y = applyAdjoint_(this,x)            
            % Reimplemented from parent class :class:`LinOp`.
            
            y = zeros([prod(this.sizeout), 2]) ;
            y(:,1) = x(:) ;
            y(:,2) = -1i*x(:) ;
            y = reshape(y, this.sizein) ;
		end
        function y = applyHtH_(this,x)            
            % Reimplemented from parent class :class:`LinOp`.
            
            y = zeros([prod(this.sizeout), 2]) ;
            x = reshape(x, [prod(this.sizeout), 2]) ;
            y(:,1) = x(:,1)+1i*x(:,2) ;
            y(:,2) = -1i*y(:,1) ;
            y = reshape(y, this.sizein) ;
		end
        function y = applyHHt_(~,x)            
            % Reimplemented from parent class :class:`LinOp`.
            
            y = 2*x ;
		end
        function M = makeHHt_(this)
            % Reimplemented from parent class :class:`LinOp`.
            M = LinOpDiag(this.sizeout, 2) ;
        end
        function y = applyInverse_(this,x)            
            % Reimplemented from parent class :class:`LinOp`.
            
            y = zeros([prod(this.sizeout), 2]) ;
            y(:,1) = real(x(:)) ;
            y(:,2) = imag(x(:)) ;
            y = reshape(y, this.sizein) ;
        end
        function M = makeInverse_(this)
            % Reimplemented from parent class :class:`LinOp`.
            
            M = OpC2R(this.sizeout) ;
        end
    end
end
