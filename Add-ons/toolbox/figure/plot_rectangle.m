%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to plot a list of rectangles
% 
% Created: 04/27/2018 (mm/dd/yyyy)
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
% Modified: 08/29/2019 (mm/dd/yyyy) If pix is empty, the figure is
% considered to be at scale.
% Author: Anthony Berdeu (Centre de Recherche en Astrophysique de Lyon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #fig# figure in which to plot the rectangles
%
% #list_center# list of the centers of the rectangles
%
% #side = [l, L]# lengths of the sides of the rectangles
%
% #theta# rotation angle of the rectangles
%
% #pix# pixel chracteristic
%   pix.nb_x -> number of the pixel on the x-axis
%   pix.nb_y -> number of the pixel on the y-axis
%   pix.dx -> size of the pixel on the x-axis
%   pix.dy -> size of the pixel on the y-axis
%
% #flag_color# color of the hexagons
%
% #LineWidth# Line width
%
%%%%%%%%
% Ouput
%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_rectangle(fig, list_center, side, theta, pix, flag_color, ...
    LineWidth)

    %% Options
    if nargin < 7 || isempty(LineWidth)
        LineWidth = 1 ;
    end
    
    %% Positions of the centers and perimeters on the picture
    nb_rect = size(list_center,1) ;
    if isempty(pix)
        x_c = list_center(:,2) ;
        y_c = list_center(:,1) ;

        x_side = ...
            [-side(1)/2, -side(1)/2, side(1)/2, side(1)/2 -side(1)/2, ...
            -side(1)/2] ;
        y_side = ...
            [-side(2)/2, side(2)/2, side(2)/2, -side(2)/2 -side(2)/2, ...
            -side(2)/2] ;

        % Rotation
        [x_side, y_side] = rot_2D(theta, x_side, y_side) ;

        x_side = meshgrid(x_side, 1:nb_rect)+x_c ;
        y_side = meshgrid(y_side, 1:nb_rect)+y_c ;
    else
        x_c = list_center(:,2)/pix.dx+floor(pix.nb_x/2)+1 ;
        y_c = list_center(:,1)/pix.dy+floor(pix.nb_y/2)+1 ;

        x_side = ...
            [-side(1)/2, -side(1)/2, side(1)/2, side(1)/2 -side(1)/2, ...
            -side(1)/2] ;
        y_side = ...
            [-side(2)/2, side(2)/2, side(2)/2, -side(2)/2 -side(2)/2, ...
            -side(2)/2] ;

        % Rotation
        [x_side, y_side] = rot_2D(theta, x_side, y_side) ;

        x_side = meshgrid(x_side,1:nb_rect)/pix.dx+x_c ;
        y_side = meshgrid(y_side,1:nb_rect)/pix.dy+y_c ;
    end
    
    
    
    %% Plot the hexagons
    figure(fig) ;
    hold on
    plot(x_c, y_c, '.', 'Color', flag_color, 'LineWidth', LineWidth) ;
    for h = 1:nb_rect
        plot(x_side(h,:),y_side(h,:), 'Color', flag_color, ...
            'LineWidth', LineWidth) ;
    end
    hold off
end