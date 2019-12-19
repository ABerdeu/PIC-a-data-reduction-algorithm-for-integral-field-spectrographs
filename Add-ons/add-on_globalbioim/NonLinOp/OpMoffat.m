classdef OpMoffat < Map
    % OpMoffat: Operator to simulate a Moffat pattern
    % 
    % :param y: list of the coordinates on the y-axis on which the Moffat
    % pattern is simulated
    % 
    % :param x: list of the coordinates on the x-axis on which the Moffat
    % pattern is simulated
    % 
    % :param y_c: coordinate of the center of the Moffat pattern on the
    % y-axis (if empty, it becomes a variable to give as an input)
    % 
    % :param x_c: coordinate of the center of the Moffat pattern on the
    % x-axis (if empty, it becomes a variable to give as an input)
    % 
    % :param alpha: elongation of the Moffat pattern along the radial
    % direction  (if empty, it becomes a variable to give as an input)
    % 
    % :param beta: power of the attenuation of the Moffat pattern (if
    % empty, it becomes a variable to give as an input)
    % 
    % :param amp: amplitude of the Moffat pattern (if empty, it becomes a
    % variable to give as an input)
    % 
    % :param offset: offset of the Moffat pattern (if empty, it becomes a
    % variable to give as an input)
    % 
    % :param flag_norm: is the Moffat pattern normalized? (used to
    % uncouple the variables 'amplitude', 'alpha' and 'beta' in the
    % optimization process) (default: false)
    % 
    % **Example** Moff = OpMoffat(y, x, y_c, x_c, alpha, beta, amp, offset)
    %
    % **Example** Moff = OpMoffat(y, x, y_c, x_c, alpha, beta, amp, offset, true)
    %
    % See also :class:`Map`
  
    %%    Copyright (C) 2018
    %     Created: 08/31/2018 (mm/dd/yyyy)
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

    % Protected Set and public Read properties
    properties (SetAccess = protected,GetAccess = public)
        y;          % list of the coordinates on the y-axis on which the
            % Moffat pattern is simulated
        x;          % list of the coordinates on the x-axis on which the
            % Moffat pattern is simulated
        list_input; % List of the parameters which are inputs
        list_par;   % List of the parameters which are fixed
        flag_norm;  % Is the Moffat pattern normalized?
    end
    
    %% Constructor
    methods
        function this = OpMoffat(varargin)
            % varargin = 
            %   {y, x, y_c, x_c, alpha, beta, amp, offset}
            this.name = 'OpMoffat' ;
            this.isDifferentiable=true ;
            this.isInvertible=false ;
            
            switch length(varargin)
                case 8
                    this.flag_norm = false ;
                case 9
                    this.flag_norm = varargin{9} ;
                otherwise
                    error(['The number of input arguments is not', ...
                        ' correct...']) ;
            end
            
            % Gaussian position
            [this.x, this.y] = meshgrid(varargin{2}, varargin{1}) ;
            this.sizeout = size(this.x) ;
            this.list_input = [] ;
            this.list_par = zeros(6,1) ;
            
            % Defining inputs and parameters
            for p = 1:6
                par_p = varargin{p+2} ;
                if isempty(par_p)
                    this.list_input = [this.list_input ; p] ;
                else
                    this.list_par(p) = par_p ;
                end
            end
            this.sizein = [length(this.list_input), 1] ;
		end
    end
	
    %% Core Methods containing implementations (Protected)
	methods (Access = protected)
        function y = apply_(this,x)
            % Reimplemented from parent class :class:`Map`.
            
            % Initialization of the parameters
            % y_c, x_c, alpha, beta, amp, offset
            par = this.list_par ;
            par(this.list_input) = x ;
            
            % Moffat simulation
            y = par(5)*( ...
                1 + ((this.x-par(2)).^2 + (this.y-par(1)).^2)/par(3).^2 ...
                ).^-par(4) ;
            
            % Is the pattern normalized?
            if this.flag_norm
                y = y*(par(4)-1)/(pi*par(3).^2) ;
            end

            % Offset
            y = y + par(6) ;
        end
        
        function x = applyJacobianT_(this,y,v)
            % Reimplemented from parent class :class:`Map`.
            
            % Initialization of the parameters
            % y_c, x_c, alpha, beta, amp, offset
            par = this.list_par ;
            par(this.list_input) = v ;
            dx = this.x(:)-par(2) ;
            dy = this.y(:)-par(1) ;
            A = 1 + (dx.^2 + dy.^2)/par(3).^2 ;
            % Is the pattern normalized?
            if this.flag_norm
                c = (par(4)-1)/(pi*par(3).^2) ;
            else
                c = 1 ;
            end
            x = zeros(this.sizein) ;
            
            for p = 1:length(this.list_input)
                switch this.list_input(p)
                    case 1 % y_c
                        x(p) = sum( ...
                            2*par(5)*par(4)/par(3)^2.*c*dy.* ...
                            A.^(-par(4)-1).*y(:)) ;
                    case 2 % x_c
                        x(p) = sum( ...
                            2*par(5)*par(4)/par(3)^2.*c*dx.* ...
                            A.^(-par(4)-1).*y(:)) ;
                    case 3 % alpha
                        if this.flag_norm
                            x(p) = sum( ...
                                2*par(5)*c/par(3).*( ...
                                par(4)/par(3).^2.*(dx.^2 + dy.^2)-A).* ...
                                A.^(-par(4)-1).*y(:)) ;
                        else
                            x(p) = sum( ...
                                2*par(5)*par(4)/par(3).^3.* ...
                                (dx.^2 + dy.^2).*A.^(-par(4)-1).*y(:)) ;
                        end
                    case 4 % beta
                        if this.flag_norm
                            x(p) = sum( ...
                                par(5)/(pi*par(3).^2).*( ...
                                1-(par(4)-1).*log(A)).* ...
                                A.^-par(4).*y(:)) ;
                        else
                            x(p) = sum(-par(5)*log(A).*A.^-par(4).*y(:)) ;
                        end
                    case 5 % amp
                        x(p) = sum(c*A.^-par(4).*y(:)) ;
                    case 6 % offset
                        x(p) = sum(y(:)) ;
                end
            end
        end
    end
end
