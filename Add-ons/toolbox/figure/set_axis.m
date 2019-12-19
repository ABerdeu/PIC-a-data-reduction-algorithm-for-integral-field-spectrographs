%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to set the axis of a figure
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
% #ax# axis to set
%   'x' -> x-axis
%   'y' -> y-axis
%
% #LineWidth# the line width of the sticks
%
% #list_tick_pos# the list of the position of the ticks
%
% #list_tick_label# the list of the labels of the ticks
%
% #FontSize_tick# the font size of the axis ticks
%
% #axis_label# label of the axis
%
% #FontSize_label# the font size of the axis label
%
% #title_label# title of the figure
%
% #FontSize_title# the font size of the figure
%
%%%%%%%%
% Ouput
%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function set_axis(h_fig, ax, LineWidth, list_tick_pos, list_tick_label, ...
    FontSize_tick, axis_label, FontSize_label, title_label, FontSize_title)
    
    % Selecting figure
    figure(h_fig) ;

    % Getting axes
    axes = gca ;

    % Line width of the axis
    axes.LineWidth = LineWidth;
    
    % Ticks
    switch ax
        case 'x'
            set(axes, 'xtick', list_tick_pos) ;
            set(axes, 'xticklabels', list_tick_label) ;
            h = get(axes, 'XAxis') ;
            set(h, 'FontSize', FontSize_tick);
        case 'y'
            set(axes, 'ytick', list_tick_pos) ;
            set(axes, 'yticklabels', list_tick_label) ;
            h = get(axes, 'YAxis') ;
            set(h, 'FontSize', FontSize_tick) ;
        otherwise
            error([axis_label, ' is an unknown axe...']) ;
    end
%     set(axes,'FontSize', FontSize_tick);
    
    % Label
    if ~isempty(axis_label)
        switch ax
            case 'x'
                xlabel(axis_label, 'FontSize', FontSize_label) ;
            case 'y'
                ylabel(axis_label, 'FontSize', FontSize_label) ;
            otherwise
                error([axis_label, ' is an unknown axe...']) ;
        end
    end
    
    % Title
    if ~isempty(title_label)
        title(title_label, 'FontSize', FontSize_title) ;
    end
    
    % Latex interpreter
    set(axes,'TickLabelInterpreter','latex') ;
    set(findall(gcf,'-property','Interpreter'), 'Interpreter', ...
        'Latex') ;

end