%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to plot a list of segments
% 
% Created: 04/03/2019 (mm/dd/yyyy)
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #fig# figure in which to plot the rectangles
%
% #list_corner# list of the corners of the polygone
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
function plot_polygone(fig, list_corner, flag_color, LineWidth)

    %% Options
    if nargin < 4 || isempty(LineWidth)
        LineWidth = 1 ;
    end
    
    %% Positions of the corners
    x_side = list_corner(:,1) ;
    x_side = [x_side; x_side(1:2)] ;
    y_side = list_corner(:,2) ;
    y_side = [y_side; y_side(1:2)] ;
    
    %% Plot the hexagons
    figure(fig) ;
    hold on
    plot(x_side, y_side, '-', 'Color', flag_color, ...
        'LineWidth', LineWidth) ;
    hold off
end