classdef OpPSF_Hex < Map
    % OpPSF_Hex: Operator to simulate the incoherent point spread function
    % of a hexagonal aperture
    % 
    % :param y: list of the coordinates on the y-axis on which the pattern
    % is simulated
    % 
    % :param x: list of the coordinates on the x-axis on which the pattern
    % is simulated
    % 
    % :param y_c: coordinate of the center of the pattern on the
    % y-axis (if empty, it becomes a variable to give as an input)
    % 
    % :param x_c: coordinate of the center of the pattern on the
    % x-axis (if empty, it becomes a variable to give as an input)
    % 
    % :param theta: orientation (in degree) of the pattern
    % (if empty, it becomes a variable to give as an input)
    % 
    % :param sig: elongation of the pattern
    %
    % :param amp: amplitude of the pattern (if empty, it becomes a
    % variable to give as an input)
    % 
    % :param offset: offset of the pattern (if empty, it becomes a
    % variable to give as an input)
    % 
    % :param eps: the minimal value to avoid a division by 0 
    % (default 1e-12)
    % 
    % **Example** Gauss = OpPSF_Hex(y, x, y_c, x_c, theta, sig, amp, offset)
    % 
    % **Example** Gauss = OpPSF_Hex(y, x, y_c, x_c, theta, sig, amp, offset, eps)
    % 
    % See also :class:`Map`
  
    %%    Copyright (C) 2018
    %     Created: 02/22/2019 (mm/dd/yyyy)
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
            % pattern is simulated
        x;          % list of the coordinates on the x-axis on which the
            % pattern is simulated
        list_input; % List of the parameters which are inputs
        list_par;   % List of the parameters which are fixed
        eps;        % The value of epsilon to avoid a division by zero
    end
    
    %% Constructor
    methods
        function this = OpPSF_Hex(varargin)
            % varargin = 
            %   {y, x, y_c, x_c, theta, sig, amp, offset, (eps)}
            this.name = 'OpPSF_Hex' ;
            this.isDifferentiable=true ;
            this.isInvertible=false ;
            
            % Checking number of output
            if length(varargin) == 9
                this.eps = varargin{9} ;
            else
                this.eps = 1e-12 ;
            end
            
            % Positions
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
            % y_c, x_c, theta, sig, amp, offset
            par = this.list_par ;
            par(this.list_input) = x ;
            
            % Sine and cosine (in degree)
            c = cosd(par(3)) ;
            s = sind(par(3)) ;
            
            % X and Y
            X = 2*pi/par(4).*(+c*(this.x-par(2))+s*(this.y-par(1))) ;
            Y = 2*pi/par(4).*(-s*(this.x-par(2))+c*(this.y-par(1))) ;
            
            % u_pm and v_pm
%             u_p = 1/3^0.5*(X/3+Y) ;
%             u_m = 1/3^0.5*(X/3-Y) ;
%             v_p = X+2*Y/3^0.5 ;
%             v_m = X-2*Y/3^0.5 ;
%             

            u_p = 1/3*(X+3^0.5*Y) ;
            u_m = 1/3*(X-3^0.5*Y) ;
            v_p = X+Y/3^0.5 ;
            v_m = X-Y/3^0.5 ;
            
            v = v_p + v_m ;
            
            % Preveting any division by zero...
            u_p(u_p==0) = this.eps ;
            sel_u_p = abs(u_p)<this.eps ;
            u_p(sel_u_p) = this.eps*sign(u_p(sel_u_p)) ;
            
            u_m(u_m==0) = this.eps ;
            sel_u_m = abs(u_m)<this.eps ;
            u_m(sel_u_m) = this.eps*sign(u_m(sel_u_m)) ;
            
            v(v==0) = this.eps ;
            sel_v = abs(v)<this.eps ;
            v(sel_v) = this.eps*sign(v(sel_v)) ;
            
            % Fourier transform of a hexagon
            H = 2./v.*( ...
                sin(v_p).*sin(u_m)./u_m + ...
                sin(v_m).*sin(u_p)./u_p) ;
            
            % Normalization, amplitude and offset
            y = par(5)*2/(3^0.5*par(4).^2)*H.^2+par(6) ;

        end
        
        function x = applyJacobianT_(this,y,v)
            % Reimplemented from parent class :class:`Map`.
            
            % Initialization of the parameters
            % y_c, x_c, theta, sig, amp, offset
            par = this.list_par ;
            par(this.list_input) = v ;

            % Sine and cosine (in degree)
            c = cosd(par(3)) ;
            s = sind(par(3)) ;
            
            % X and Y (vector shape)
            X = 2*pi/par(4).*(+c*(this.x-par(2))+s*(this.y-par(1))) ;
            Y = 2*pi/par(4).*(-s*(this.x-par(2))+c*(this.y-par(1))) ;
            X = X(:) ;
            Y = Y(:) ;
            
            % u_pm and v_pm
            u_p = 1/3*(X+3^0.5*Y) ;
            u_m = 1/3*(X-3^0.5*Y) ;
            v_p = X+Y/3^0.5 ;
            v_m = X-Y/3^0.5 ;
            
            v = v_p + v_m ;
            
            % Preveting any division by zero...
            u_p(u_p==0) = this.eps ;
            sel_u_p = abs(u_p)<this.eps ;
            u_p(sel_u_p) = this.eps*sign(u_p(sel_u_p)) ;
            
            u_m(u_m==0) = this.eps ;
            sel_u_m = abs(u_m)<this.eps ;
            u_m(sel_u_m) = this.eps*sign(u_m(sel_u_m)) ;
            
            v(v==0) = this.eps ;
            sel_v = abs(v)<this.eps ;
            v(sel_v) = this.eps*sign(v(sel_v)) ;
            
            % Fourier transform of a hexagon
            H = 2./v.*( ...
                sin(v_p).*sin(u_m)./u_m + ...
                sin(v_m).*sin(u_p)./u_p) ;
            
            % Main derivative
            dH_du_m = 2./v.*sin(v_p)./u_m.^2.*(u_m.*cos(u_m)-sin(u_m)) ;
            dH_du_p = 2./v.*sin(v_m)./u_p.^2.*(u_p.*cos(u_p)-sin(u_p)) ;
            dH_dv_m = 1./v.*(2*sin(u_p)./u_p.*cos(v_m)-H) ;
            dH_dv_p = 1./v.*(2*sin(u_m)./u_m.*cos(v_p)-H) ;
            
            % Norm factor
            A = par(5)*4/(3^0.5*par(4)^2) ;
            
            x = zeros(this.sizein) ;

            % Loop on the variables
            for p = 1:length(this.list_input)
                switch this.list_input(p)
                    case 1 % y_c
                        dH_dy = ( ...
                                2*pi/(3*par(4))*(-s+3^0.5*c).*dH_du_m + ...
                                2*pi/(3*par(4))*(-s-3^0.5*c).*dH_du_p + ...
                                2*pi/par(4)*(-s+1/3^0.5*c).*dH_dv_m + ...
                                2*pi/par(4)*(-s-1/3^0.5*c).*dH_dv_p ...
                            ) ;
                        x(p) = sum(A*dH_dy.*H.*y(:)) ;
                    case 2 % x_c
                        dH_dx = ( ...
                                2*pi/(3*par(4))*(-c-3^0.5*s).*dH_du_m + ...
                                2*pi/(3*par(4))*(-c+3^0.5*s).*dH_du_p + ...
                                2*pi/par(4)*(-c-1/3^0.5*s).*dH_dv_m + ...
                                2*pi/par(4)*(-c+1/3^0.5*s).*dH_dv_p ...
                            ) ;
                        x(p) = sum(A*dH_dx.*H.*y(:)) ;
                    case 3 % theta
                        dH_dtheta = ( ...
                                1/3*pi/180.*(Y+3^0.5*X).*dH_du_m + ...
                                1/3*pi/180.*(Y-3^0.5*X).*dH_du_p + ...
                                pi/180.*(Y+1/3^0.5*X).*dH_dv_m + ...
                                pi/180.*(Y-1/3^0.5*X).*dH_dv_p ...
                            ) ;
                        x(p) = sum(A*dH_dtheta.*H.*y(:)) ;
                    case 4 % sigma
                        dH_sigma = ( ...
                                -1/par(4).*u_m.*dH_du_m + ...
                                -1/par(4).*u_p.*dH_du_p + ...
                                -1/par(4).*v_m.*dH_dv_m + ...
                                -1/par(4).*v_p.*dH_dv_p ...
                            ) ;
                        x(p) = sum(A*(dH_sigma-1/par(4)*H).*H.*y(:)) ;
                    case 5 % amp
                        x(p) = sum(2/(3^0.5*par(4)^2)*H.^2.*y(:)) ;
                    case 6 % offset
                        x(p) = sum(y(:)) ;
                end
            end
        end
    end
end
