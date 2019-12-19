%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to simulation a guassian pattern according to its position and 
% parameters
% 
% Created: 05/15/2018 (mm/dd/yyyy)
% Modified: 05/17/2018 (mm/dd/yyyy) More possibilities for the gaussian
% simulation
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #[list_y, list_x]# list of the position of the gaussian
%   pattern in the figure
%
% #gauss_pos = [y_c, x_c]# center of the gaussian pattern
%
% #par# parameters of the gaussian pattern
%   -> sig_x: gaussian extension on the x-axes
%   -> sig_y: gaussian extension on the y-axes
%   -> amp: amplitude of the gaussian patterns (default 1)
%
% #theta# orientation of the gaussian pattern (deg)
%
%%%%%%%%
% Ouput
%%%%%%%%
% #gauss_pic# the picture of the simulated gaussian pattern
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gauss_pic = get_gauss_simulation(list_y, list_x, gauss_pos, ...
    par, theta)

    %% Initialization
    sig_x = par(1) ;
    sig_y = par(2) ;
    if length(par) == 3
        amp = par(3) ;
    else
        amp = 1 ;
    end
    if nargin < 4
        theta = 0 ;
    end
    
    %% Rotation on the frame
    [list_xx, list_yy] = meshgrid(list_x, list_y) ;
    list_xx = list_xx-gauss_pos(2) ;
    list_yy = list_yy-gauss_pos(1) ;
    [list_xx, list_yy] = rot_2D(-theta, list_xx, list_yy) ;

    % Computation of the elementary PSF
    gauss_pic = amp*exp(-0.5*((list_xx/sig_x).^2 + (list_yy/sig_y).^2)) ;
end