classdef OpGauss < Map
    % OpGauss: Operator to simulate a gaussian pattern
    % 
    % :param y: list of the coordinates on the y-axis on which the gaussian
    % pattern is simulated
    % 
    % :param x: list of the coordinates on the x-axis on which the gaussian
    % pattern is simulated
    % 
    % :param y_c: coordinate of the center of the gaussian pattern on the
    % y-axis (if empty, it becomes a variable to give as an input)
    % 
    % :param x_c: coordinate of the center of the gaussian pattern on the
    % x-axis (if empty, it becomes a variable to give as an input)
    % 
    % :param theta: orientation (in degree) of the gaussian pattern
    % (if empty, it becomes a variable to give as an input)
    % This field is not used in the case of axisymmetric pattern.
    % 
    % :param sig_par: elongation of the gaussian pattern along the
    % direction given by theta (if empty, it becomes a variable to give as
    % an input)
    % This field is not used in the case of axisymmetric pattern. It is
    % replaced by sigma, the radial elongation of the gaussian pattern.
    % 
    % :param sig_perp: elongation of the gaussian pattern along the
    % perpendicular direction given by theta (if empty, it becomes a
    % variable to give as an input)
    % This field is not used in the case of axisymmetric pattern. It is
    % replaced by sigma, the radial elongation of the gaussian pattern.
    % 
    % :param amp: amplitude of the gaussian pattern (if empty, it becomes a
    % variable to give as an input)
    % 
    % :param offset: offset of the gaussian pattern (if empty, it becomes a
    % variable to give as an input)
    % 
    % :param flag_norm: is the gaussian pattern normalized? (used to
    % uncouple the variables 'amplitude' and 'sigma' in the optimization
    % process) (default: false)
    % 
    % **Example** Gauss = OpGauss(y, x, y_c, x_c, theta, sig_par, sig_perp, amp, offset)
    % 
    % **Example** Gauss = OpGauss(y, x, y_c, x_c, theta, sig_par, sig_perp, amp, offset, flag_norm)
    % 
    % **Example** Gauss = OpGauss(y, x, y_c, x_c, sig, amp, offset)
    %
    % **Example** Gauss = OpGauss(y, x, y_c, x_c, sig, amp, offset, flag_norm)
    %
    % See also :class:`Map`
  
    %%    Copyright (C) 2018
    %     Created: 07/26/2018 (mm/dd/yyyy)
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
            % gaussian pattern is simulated
        x;          % list of the coordinates on the x-axis on which the
            % gaussian pattern is simulated
        flag_sym ;  % Is the pattern axisymmetric?
        list_input; % List of the parameters which are inputs
        list_par;   % List of the parameters which are fixed
        flag_norm;  % Is the gaussian pattern normalized?
    end
    
    %% Constructor
    methods
        function this = OpGauss(varargin)
            % varargin = 
            %   {y, x, y_c, x_c, theta, sig_par, sig_perp, amp, offset}
            this.name = 'OpGauss' ;
            this.isDifferentiable=true ;
            this.isInvertible=false ;
            
            switch length(varargin)
                case 7
                    this.flag_sym = true ;
                    this.flag_norm = false ;
                case 8
                    this.flag_sym = true ;
                    this.flag_norm = varargin{8} ;
                case 9
                    this.flag_sym = false ;
                    this.flag_norm = false ;
                case 10
                    this.flag_sym = false ;
                    this.flag_norm = varargin{10} ;
                otherwise
                    error(['The number of input arguments is not', ...
                        ' correct...']) ;
            end
            
            % Gaussian position (the problem is separable on x and y for 
            % symmetric patterns
            % y*x'
            this.y = varargin{1}(:) ; % Vector shape
            this.x = varargin{2}(:) ; % Vector shape
            this.sizeout = [length(this.y), length(this.x)] ;
            if ~this.flag_sym
                [this.x, this.y] = meshgrid(this.x, this.y) ;
            end
            this.list_input = [] ;
            this.list_par = zeros((length(varargin)-2),1) ;
            
            % Defining inputs and parameters
            for p = 1:(length(varargin)-2)
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
            
            if this.flag_sym
                % Initialization of the parameters
                % y_c, x_c, sig, amp, offset
                par = this.list_par ;
                par(this.list_input) = x ;

                % Gaussian simulation
                exp_y = exp(-0.5*(this.y-par(1)).^2 / par(3).^2) ;
                exp_x = exp(-0.5*(this.x-par(2)).^2 / par(3).^2) ;
                y = par(4)*exp_y.*exp_x' ;

                
                % Is the pattern normalized?
                if this.flag_norm
                    y = y/(2*pi*par(3).^2) ;
                end
                
                % Offset
                y = y + par(5) ;
            else
                % Initialization of the parameters
                % y_c, x_c, theta, sig_par, sig_perp, amp, offset
                par = this.list_par ;
                par(this.list_input) = x ;
                c = cosd(par(3)) ;
                s = sind(par(3)) ;

                % Gaussian simulation
                y = par(6)*exp(-0.5*( ...
                    ((+c*(this.x-par(2))+s*(this.y-par(1)))/par(4)).^2 ...
                    + ...
                    ((-s*(this.x-par(2))+c*(this.y-par(1)))/par(5)).^2 ...
                    )) ;
                
                % Is the pattern normalized?
                if this.flag_norm
                    y = y/(2*pi*par(4)*par(5)) ;
                end
                
                % Offset
                y = y + par(7) ;
            end
        end
        
        function x = applyJacobianT_(this,y,v)
            % Reimplemented from parent class :class:`Map`.
            
            if this.flag_sym
                % Initialization of the parameters
                % y_c, x_c, sig, amp, offset
                par = this.list_par ;
                par(this.list_input) = v ;
                dx = this.x(:)-par(2) ;
                dy = this.y(:)-par(1) ;
                E_x = exp(-0.5*dx.^2/par(3).^2) ;
                E_y = exp(-0.5*dy.^2/par(3).^2) ;
                E = E_y.*E_x' ;
                
                % Is the pattern normalized?
                if this.flag_norm
                    E = E/(2*pi*par(3).^2) ;
                end
                x = zeros(this.sizein) ;
                
                % Loop on the variables
                for p = 1:length(this.list_input)
                    switch this.list_input(p)
                        case 1 % y_c
                            E_p = dy.*E ;
                            x(p) = sum(E_p(:)./par(3)^2.*par(4).*y(:)) ;
                        case 2 % x_c
                            E_p = E.*dx' ;
                            x(p) = sum(E_p(:)./par(3)^2.*par(4).*y(:)) ;
                        case 3 % sig
                            % Is the pattern normalized?
                            if this.flag_norm
                                E_p = (dy.^2+dx'.^2) ;
                                x(p) = sum((E_p(:)./par(3)^3 - ...
                                    2/par(3)).*par(4).*E(:).*y(:)) ;
                            else
                                E_p = (dy.^2+dx'.^2) ;
                                x(p) = sum(E_p(:)./par(3)^3.* ...
                                    par(4).*E(:).*y(:)) ;
                            end
                        case 4 % amp
                            x(p) = sum(E(:).*y(:)) ;
                        case 5 % offset
                            x(p) = sum(y(:)) ;
                    end
                end
            else
                % Initialization of the parameters
                % y_c, x_c, theta, sig_par, sig_perp, amp, offset
                par = this.list_par ;
                par(this.list_input) = v ;
                c = cosd(par(3)) ;
                s = sind(par(3)) ;
                inv_sig_par = 1/par(4) ;
                inv_sig_perp = 1/par(5) ;
                delta_inv_sig2 = inv_sig_par^2-inv_sig_perp^2 ;
                dx = this.x(:)-par(2) ;
                dy = this.y(:)-par(1) ;
                E = exp(-0.5*( ...
                    ((+c*dx + s*dy)/par(4)).^2 + ...
                    ((-s*dx + c*dy)/par(5)).^2 ...
                    )) ;
                % Is the pattern normalized?
                if this.flag_norm
                    E = E/(2*pi*par(4)*par(5)) ;
                end
                x = zeros(this.sizein) ;

                % Loop on the variables
                for p = 1:length(this.list_input)
                    switch this.list_input(p)
                        case 1 % y_c
                            x(p) = sum( ...
                                ((inv_sig_perp^2+delta_inv_sig2*s^2).* ...
                                    dy + ...
                                0.5*dx*sin(pi/90*par(3)).* ...
                                    delta_inv_sig2 ...
                                )*par(6).*E.*y(:)) ;
                        case 2 % x_c
                            x(p) = sum( ...
                                ((inv_sig_perp^2+delta_inv_sig2*c^2).* ...
                                    dx + ...
                                0.5*dy*sin(pi/90*par(3)).* ...
                                    delta_inv_sig2 ...
                                )*par(6).*E.*y(:)) ;
                        case 3 % theta
                            x(p) = pi/360*delta_inv_sig2*sum(...
                                ( ...
                                sin(pi/90*par(3)) * (dx.^2 + dy.^2) ...
                                - 2*cos(pi/90*par(3)) *dx.*dy ...
                                )*par(6).*E.*y(:)) ;
                        case 4 % sig_par
                            % Is the pattern normalized?
                            if this.flag_norm
                                x(p) = sum( ...
                                    (inv_sig_par^3*(+c*dx + s*dy).^2 - ...
                                    inv_sig_par).*par(6).*E(:).*y(:)) ;
                            else
                                x(p) = sum( ...
                                    inv_sig_par^3*(+c*dx + s*dy).^2.* ...
                                    par(6).*E(:).*y(:)) ;
                            end
                        case 5 % sig_perp
                            if this.flag_norm
                                x(p) = sum( ...
                                    (inv_sig_perp^3*(-s*dx + c*dy).^2 - ...
                                    inv_sig_perp).*par(6).*E(:).*y(:)) ;
                            else
                                x(p) = sum( ...
                                    inv_sig_perp^3*(-s*dx + c*dy).^2* ...
                                    par(6).*E(:).*y(:)) ;
                            end
                        case 6 % amp
                            x(p) = sum(E(:).*y(:)) ;
                        case 7 % offset
                            x(p) = sum(y(:)) ;
                    end
                end
            end
        end
    end
end
