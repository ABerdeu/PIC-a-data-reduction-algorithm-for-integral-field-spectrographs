%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get the linear convolution operator of a hexagonal aperture
% 
% Created: 04/04/2018 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #pix# pixel chracteristic
%   pix.dx -> size of the pixel on the x-axis
%   pix.dy -> size of the pixel on the y-axis
%   pix.nb_x_pad -> number of the pixel on the x-axis in the padded domain
%   pix.nb_y_pad -> number of the pixel on the y-axis in the padded domain
%
% #side# side of of the hexagon compared to the cartesian grid
% of the lenslet grid
%
% #theta# orientation angle of the hexagon compared to the cartesian grid
% of the lenslet grid
%
% #nb_l# number of wavelength one which to perform the convolution
%
% #flag_space#
%   SD -> compute the convolution kernel in the spatial domain
%   FD -> compte the convolution kernel in the Fourier domain
%
% #const#
%   flag_space = SD -> factor of scaling to compute the apodization on the
%   edges of the hexagon
%   flag_space = FD -> numerical threshold to use Taylor expansion of the
%       formula (if 0, the numerical computation is supposed to be valid)
%
%%%%%%%%
% Ouput
%%%%%%%%
% #Conv_hexa# the convolution linear operator for the hexagonal aperture
%
% #FT_hexagon# the Fourier transform of the hexagon
%
% #flag_norm# Normalization of the convolution? (default 0)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Conv_hexa, FT_hexagon] = get_conv_hexagon(pix, side, theta, ...
    nb_l, flag_space, const, flag_norm)

    if nargin < 7
        flag_norm = false ;
    end

    switch flag_space
        case 'SD'
            %% Initialization
            % Initialization of the impulse response
            IR = zeros([pix.nb_y_pad, pix.nb_x_pad]) ;

            % Parameters of the scaled hexagon
            const = floor(const)+1 ;
            pix_scaled.nb_y = (floor(2*(side+2*pix.dy)*const/pix.dy)+1) ;
            pix_scaled.nb_x = (floor(2*(side+2*pix.dx)*const/pix.dx)+1) ;
            pix_scaled.dy = pix.dy/const ;
            pix_scaled.dx = pix.dx/const ;

            % Coordinates in the scaled spatial domain
            [x_2D, y_2D] = meshgrid( ...
                get_Fourier_vector(pix_scaled.nb_x, pix_scaled.dx), ...
                get_Fourier_vector(pix_scaled.nb_y, pix_scaled.dy)) ;


            %% Computation of the impulse response in the scaled frame
            IR_scaled = true(pix_scaled.nb_y, pix_scaled.nb_x) ;
            % Threshold for the hexagon edges
            th_y = 3^0.5*side/2 ;

            % Loop on the three parts of the hexagon
            for theta_hex = 0:60:120
                % Rotation of the coordinates
                [~, r_y] = rot_2D(theta_hex-theta, x_2D, y_2D) ;

                % Mask on the part of the hexagon
                IR_scaled = IR_scaled & (abs(r_y)<th_y) ;
            end

            %% Rescaling
            list_ind = get_Fourier_vector(const, 1) ;
            size_y = floor(side/pix.dy)+1 ;
            size_x = floor(side/pix.dx)+1 ;

            % Loop on the positions
            for iy = -size_y:size_y
                for jx = -size_x:size_x
                    % Position in the scaled matrix
                    ij_scaled = pos2ind([iy*pix.dy, jx*pix.dx], ...
                        pix_scaled) ;

                    % Extraction of the scaled area
                    c = IR_scaled(ij_scaled(1)+list_ind, ...
                        ij_scaled(2)+list_ind) ;

                    % Update of the impulse response
                    IR(iy+floor(pix.nb_y_pad/2)+1, ...
                        jx+floor(pix.nb_x_pad/2)+1) = mean(c(:)) ;
                end
            end        

            %% Operator declaration
            if flag_norm
                IR = IR/sum(IR(:)) ;
            end
            
            FT_hexagon = fftshift(fft2(ifftshift(IR))) ;
            IR = bsxfun(@times, ifftshift(IR), ...
                ones(pix.nb_y_pad, pix.nb_x_pad, nb_l));
            Conv_hexa = LinOpConv('PSF', IR, true, 1:2) ;

        case 'FD'
            %% Initialization
            % Coordinates in the Fourier domain
            [u_2D, v_2D] = meshgrid( ...
                get_Fourier_vector(pix.nb_x_pad, ...
                    1/(pix.nb_x_pad*pix.dx)), ...
                get_Fourier_vector(pix.nb_y_pad, ...
                    1/(pix.nb_y_pad*pix.dy))) ;

            % Rotation
            [u_2D, v_2D] = rot_2D(-theta, u_2D, v_2D) ;

            %% Computation of the transfert function in the Fourier domain
            side_hex = 3^0.5*side/2 ;
            b_0 = 4*pi*side_hex*u_2D/3^0.5 ;
            b_p = 2*pi*side_hex*(u_2D/3^0.5 + v_2D) ;
            b_m = 2*pi*side_hex*(u_2D/3^0.5 - v_2D) ;

            % Global formula
            TF = 4*side_hex^2* ...
                (b_m.*cos(b_m) - b_0.*cos(b_0) + b_p.*cos(b_p)) ...
                ./(3^0.5*b_m.*b_0.*b_p);

            % Taking care of the 0-values
            epsilon = pi*side_hex*v_2D ;
            gamma = 4*epsilon ;

            % u=0
            if const>0
                ind = abs(u_2D) < const ;
                TF( ind ) = 2*side_hex^2/3^0.5.* ...
                    (2*cos(epsilon(ind))+sin(epsilon(ind))./ ...
                    epsilon(ind)).*sin(epsilon(ind))./epsilon(ind) ;
            else
                ind = u_2D==0 ;
                TF( ind ) = side_hex^2*( ...
                    gamma(ind).*cos(epsilon(ind)) + ...
                    2*sin(epsilon(ind))) ...
                    .*sin(epsilon(ind))./(3^0.5*epsilon(ind).^2) ;
            end

            % u=+3^0.5*v
            if const>0
                ind = abs(u_2D-3^0.5*v_2D) < const ;
                TF( ind ) = 4*side_hex^2/3^0.5.* ...
                    (gamma(ind).^2/2+sin(epsilon(ind))./epsilon(ind)) ;
            else
                ind = u_2D==3^0.5*v_2D ;
                TF(ind) = 4*side_hex^2* ...
                    (1-cos(gamma(ind)) + ...
                    gamma(ind)* ...
                    sin(gamma(ind))) ...
                    ./(3^0.5*gamma(ind).^2) ;
            end

            % u=-3^0.5*v
            if const>0
                ind = abs(u_2D+3^0.5*v_2D) < const ;
                TF( ind ) = 4*side_hex^2/3^0.5.* ...
                    (gamma(ind).^2/2+sin(epsilon(ind))./epsilon(ind)) ;
            else
                ind = u_2D==-3^0.5*v_2D ;
                TF(ind) = 4*side_hex^2* ...
                    (1-cos(gamma(ind)) + ...
                    gamma(ind)* ...
                    sin(gamma(ind))) ...
                    ./(3^0.5*gamma(ind).^2) ;
            end

            % u=0 and v=0
            if const>0
                TF( abs(u_2D)<const & abs(v_2D)<const ) = ...
                    2*3^0.5*side_hex^2 ;
            else
                TF( (u_2D==0) & (v_2D==0) ) = 2*3^0.5*side_hex^2 ;
            end

            %% Operator declaration
            FT_hexagon = TF ;
            TF = bsxfun(@times, ifftshift(TF), ...
                ones(pix.nb_y_pad, pix.nb_x_pad, nb_l));
            Conv_hexa = LinOpConv('MTF', TF, true, 1:2) ;
            
            if flag_norm
                norm_factor = zeros(Conv_hexa.sizein) ;
                norm_factor(1,1) = 1 ;
                norm_factor = fftshift(norm_factor) ;
                norm_factor = Conv_hexa_2D*norm_factor ;
                norm_factor = sum(norm_factor(:)) ;
                Conv_hexa = LinOpConv('MTF', TF/norm_factor, true, 1:2) ;
            end

        otherwise
            error([flag_space, ': unknown model...', ...
                'Implemented models:',  ...
                newline, ...
                '   -> SD: in the spatial domain',  ...
                newline, ...
                '   -> FD: in the Fourier domain']) ;
    end
end