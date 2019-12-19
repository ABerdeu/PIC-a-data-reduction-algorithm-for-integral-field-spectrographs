%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to set the legend of a figure
% 
% Created: 07/15/2018 (mm/dd/yyyy)
% Author: Anthony Berdeu (Centre de Recherche en Astrophysique de
% Lyon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #h_fig# the figure handler
%
% #list_legend# list of the legends
%
% #FontSize# the font size
%
% #pos_legend# position of the legend
%
% #Orientation# orientation of the legend (default: 'horizontal')
%
% #EdgeColor# color of the edge (default: 'none')
%
% #Color# color of the background (default: 'none')
%
%%%%%%%%
% Ouput
%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function set_legend(h_fig, list_legend, FontSize, pos_legend, ...
    Orientation, EdgeColor, Color)
    
    % Selecting figure
    figure(h_fig) ;
    
    % Default values?
    if isempty(Orientation)
        Orientation = 'horizontal' ;
    end
    if isempty(EdgeColor)
        EdgeColor = 'none' ;
    end
    if isempty(Color)
        Color = 'none' ;
    end
    
    % Setting the position
    axes = gca ;
    lim_x = xlim ;
    lim_y = ylim ;
    frame_pos = axes.Position(1:2) ;
    frame_size = axes.Position(3:4) ;
    
    pos_legend = [frame_pos + frame_size.* ...
        (pos_legend(1:2)-[lim_x(1), lim_y(1)])./ ...
        [lim_x(2)-lim_x(1), lim_y(2)-lim_y(1)], ...
        frame_size.*pos_legend(3:4)./ ...
        [lim_x(2)-lim_x(1), lim_y(2)-lim_y(1)]] ;
    
    
    % Building legend
    legend(list_legend, ...
        'Location', pos_legend, ...
        'FontSize', FontSize, ...
        'Orientation', Orientation, ...
        'EdgeColor', EdgeColor, ...
        'Color', Color) ;
    
    % Latex interpreter
    set(findall(gcf,'-property','Interpreter'), 'Interpreter', ...
        'Latex') ;
    pause(0.1) ;
end