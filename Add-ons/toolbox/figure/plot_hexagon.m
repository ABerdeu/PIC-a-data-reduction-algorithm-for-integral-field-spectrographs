%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to plot a list of hexagon
% 
% Created: 03/21/2018 (mm/dd/yyyy)
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
% Modified: 08/29/2019 (mm/dd/yyyy) If pix is empty, the figure is
% considered to be at scale.
% Author: Anthony Berdeu (Centre de Recherche en Astrophysique de Lyon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #fig# figure in which to plot the hexagons
%
% #list_center# list of the centers of the hexagon
%
% #theta# orientation angle of the hexagons
%
% #side# side of the hexagons
%
% #pix# pixel chracteristic (optional: if not specified, the figures is
%   considered to be at scale)
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
function plot_hexagon(fig, list_center, theta, side, pix, flag_color, ...
    LineWidth)

    %% Options
    if nargin < 7 || isempty(LineWidth)
        LineWidth = 1 ;
    end

    %% Positions of the centers and edges on the picture
    nb_hexagon = size(list_center,1) ;
    if isempty(pix)
        x_c = list_center(:,2) ;
        y_c = list_center(:,1) ;

        x_side = side.* ...
            cosd(meshgrid(0:60:360,1:nb_hexagon)+theta)+x_c ;
        y_side = side.* ...
            sind(meshgrid(0:60:360,1:nb_hexagon)+theta)+y_c ;
    else
        x_c = list_center(:,2)/pix.dx+floor(pix.nb_x/2)+1 ;
        y_c = list_center(:,1)/pix.dy+floor(pix.nb_y/2)+1 ;

        x_side = side.* ...
            cosd(meshgrid(0:60:360,1:nb_hexagon)+theta)/pix.dx+x_c ;
        y_side = side.* ...
            sind(meshgrid(0:60:360,1:nb_hexagon)+theta)/pix.dy+y_c ;
    end

    %% Plot the hexagons
    figure(fig) ;
    hold on
    plot(x_c, y_c, '.', 'Color', flag_color, ...
        'LineWidth', LineWidth) ;
    for h = 1:nb_hexagon
        plot(x_side(h,:),y_side(h,:), 'Color', flag_color, ...
            'LineWidth', LineWidth) ;
    end
    hold off
end