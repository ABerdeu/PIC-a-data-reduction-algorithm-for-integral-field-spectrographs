%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to plot a list of dots
% 
% Created: 05/23/2018 (mm/dd/yyyy)
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
% #pix# pixel chracteristic (optional: if not specified, the figures is
%   considered to be at scale)
%   pix.nb_x -> number of the pixel on the x-axis
%   pix.nb_y -> number of the pixel on the y-axis
%   pix.dx -> size of the pixel on the x-axis
%   pix.dy -> size of the pixel on the y-axis
%
% #flag_color# color of the hexagons
%
% #MarkerSize# Marker size
%
%%%%%%%%
% Ouput
%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_dot(fig, list_center, pix, flag_color, MarkerSize)

    %% Options
    if nargin < 5 || isempty(MarkerSize)
        MarkerSize = 1 ;
    end
    
    %% Positions of the dots
    if isempty(pix)
        x_c = list_center(:,2) ;
        y_c = list_center(:,1) ;
    else
        x_c = list_center(:,2)/pix.dx+floor(pix.nb_x/2)+1 ;
        y_c = list_center(:,1)/pix.dy+floor(pix.nb_y/2)+1 ;
    end
    

    %% Plot the dots
    figure(fig) ;
    hold on
    plot(x_c, y_c, '.', 'Color', flag_color, 'MarkerSize', MarkerSize) ;
    hold off
end