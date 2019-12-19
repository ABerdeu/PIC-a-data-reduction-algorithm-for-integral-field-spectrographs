%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to insert a zoom in a region of interest in a picture
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
% #ROI_corner# list of two opposite corners of the region of interest
%
% #ROI_color# color of the region of interest
%
% #ROI_LineWidth# line width of the region of interest
%
% #ROI_pos# position of the inserted region of interest (x, y)
%
% #ROI_size# size of the region of interest
%
% #LineWidth# line width of the axis
%
% #FontSize# the font size of the axes ticks
%
% #ratio_tick# ratio of the tick length compared to the main figure
%
% #list_xtick_pos# the list of the position of the x-axis ticks
%
% #list_xtick_label# the list of the labels of the x-axis ticks
%
% #list_tick_ypos# the list of the position of the y-axis ticks
%
% #list_tick_ylabel# the list of the labels of the y-axis ticks
%
%%%%%%%%
% Ouput
%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function insert_ROI(h_fig, ROI_corner, ROI_color, ROI_LineWidth, ...
    ROI_pos, ROI_size, LineWidth, FontSize, ratio_tick, ...
    list_xtick_pos, list_xtick_label, ...
    list_ytick_pos, list_ytick_label)
    
    % Selecting figure
    figure(h_fig) ;
    axes_main = gca(h_fig) ;

    % Frame the ROI
    plot_polygone(h_fig, ...
        [ ...
        ROI_corner(1,1), ROI_corner(1,2) ; ...
        ROI_corner(1,1), ROI_corner(2,2) ; ...
        ROI_corner(2,1), ROI_corner(2,2) ; ...
        ROI_corner(2,1), ROI_corner(1,2) ; ...
        ], ...
        ROI_color, ROI_LineWidth) ;
    
    % Get the absolute position of the plotted frame
    lim_x = xlim ;
    lim_y = ylim ;
    frame_pos = axes_main.Position(1:2) ;
    frame_size = axes_main.Position(3:4) ;
    
    % Inserting the ROI
    axes_ROI = axes('Position', ...
        [frame_pos + frame_size.* ...
        ROI_pos./[lim_x(2)-lim_x(1), lim_y(2)-lim_y(1)], ...
        frame_size.*ROI_size./[lim_x(2)-lim_x(1), lim_y(2)-lim_y(1)]]) ;
    box on;
    
    % Copying curves
    copyobj(findobj(h_fig, 'visible', 'on','type','line'), axes_ROI);
    axis([min(ROI_corner(:,1)), max(ROI_corner(:,1)),...
        min(ROI_corner(:,2)), max(ROI_corner(:,2))]) ;
    
    % Inserting axis
    set_axis(h_fig, 'x', LineWidth, ...
        list_xtick_pos, list_xtick_label, FontSize, ...
        [], [], ...
        [], []) ;
    set_axis(h_fig, 'y', LineWidth, ...
        list_ytick_pos, list_ytick_label, FontSize, ...
        [], [], ...
        [], []) ;
    
    % Tick length
    axes_ROI.TickLength = ratio_tick*axes_main.TickLength ;
    
    % Going back to the main figure
    axes_main = axes('OuterPosition', [0, 0, 1, 1]) ;
    axis(axes_main, 'off') ;
    axis([lim_x, lim_y]) ;
end