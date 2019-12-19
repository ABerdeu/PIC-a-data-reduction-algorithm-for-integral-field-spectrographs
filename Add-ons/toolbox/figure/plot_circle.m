%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to plot a list of circles
% 
% Created: 04/04/2018 (mm/dd/yyyy)
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
% Modified: 08/29/2019 (mm/dd/yyyy) If pix is empty, the figure is
% considered to be at scale.
% Author: Anthony Berdeu (Centre de Recherche en Astrophysique de Lyon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #fig# figure in which to plot the circles
%
% #list_center# list of the centers of the circles
%
% #radii# radii of the circles
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
function plot_circle(fig, list_center, radii, pix, flag_color, LineWidth)

    %% Options
    if nargin < 6 || isempty(LineWidth)
        LineWidth = 1 ;
    end
    
    %% Positions of the centers and perimeters on the picture
    nb_circle = size(list_center,1) ;
    if isempty(pix)
        x_c = list_center(:,2) ;
        y_c = list_center(:,1) ;

        x_side = radii.*cosd(meshgrid(0:1:360,1:nb_circle))+x_c ;
        y_side = radii.*sind(meshgrid(0:1:360,1:nb_circle))+y_c ;
    else
        x_c = list_center(:,2)/pix.dx+floor(pix.nb_x/2)+1 ;
        y_c = list_center(:,1)/pix.dy+floor(pix.nb_y/2)+1 ;

        x_side = radii.*cosd(meshgrid(0:1:360,1:nb_circle))/pix.dx+x_c ;
        y_side = radii.*sind(meshgrid(0:1:360,1:nb_circle))/pix.dy+y_c ;
    end
    

    %% Plot the circles
    figure(fig) ;
    hold on
    plot(x_c, y_c, '.', 'Color', flag_color, ...
        'LineWidth', LineWidth) ;
    for h = 1:nb_circle
        plot(x_side(h,:),y_side(h,:), 'Color', flag_color, ...
            'LineWidth', LineWidth) ;
    end
    hold off
end