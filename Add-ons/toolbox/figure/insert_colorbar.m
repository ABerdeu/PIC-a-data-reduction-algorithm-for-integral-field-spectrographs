%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to insert a color bar in a picture
% 
% Created: 07/23/2018 (mm/dd/yyyy)
% Author: Anthony Berdeu (Centre de Recherche en Astrophysique de
% Lyon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #h_fig# the figure handler
%
% #color_map# color_map of the figure
%
% #flag_pos# flag on the position of the color bar
%   't' -> top
%   'b' -> bottom
%   'r' -> right
%   'l' -> left
%
% #bar_width# width of the color bar (pixels)
%
% #bar_dist# distance of the color bar from the figure (pixels)
%
% #LineWidth# line width to frame the color bar
%
% #[val_min, val_max]# extremal value of the color bar projected on 0 and 1
%
% #scale_func# scaling function (default: identity)
%
% #struct_tick# structure with the tick properties
%   .list       -> the list of the ticks
%   .LineWidth  -> the tick line width
%   .ratio      -> the tick ratio compared to the width of the color bar
%   .color      -> the color of the ticks (can be a list of color)
%
% #struct_tickLabel# structure with the tick properties
%   .list       -> the list of the tick positions
%   .listLabel  -> the list of the tick labels
%   .FontSize   -> the tick font size
%   .Alignment  -> the tick alignment (can be a list) (default: center)
%       {'left', 'center', 'right'}
%   .theta      -> the orientation (can be a list of angle) (default: 0)
%   .dist       -> distance to the color bar (default: 0)
%
%%%%%%%%
% Ouput
%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function insert_colorbar(h_fig, color_map, flag_pos, bar_width, ...
    bar_dist, LineWidth, val_min, val_max, scale_func, ...
    struct_tick, struct_tickLabel)
    
    %% Selecting figure
    figure(h_fig) ;
    axes_main = gca(h_fig) ;
    
    %% Get position and size of the figure in pixel
    % Get the new aspect ratio data
    aspect = get(axes_main,'PlotBoxAspectRatio') ;
    % Change axes Units property
    % (this only works with non-normalized units)
    set(axes_main,'Units','pixels') ;
    % Get and update the axes width to match the plot aspect ratio.
    pos = get(axes_main,'Position') ;
    pos(3) = aspect(1)/aspect(2)*pos(4);
    set(axes_main, 'Position', pos);
    
    lim_x = xlim ;
    lim_y = ylim ;

    fig_size = pos(3:4) ;
    fig_pos = pos(1:2) ;
    
    %% Creating new axes
    nb_color = size(color_map, 1) ;
    list_pos = (0:(nb_color-1))/(nb_color-1) ;
    switch flag_pos
        case 't'
            bar_pos = fig_pos ;
            bar_pos(2) = bar_pos(2)+fig_size(2)+bar_dist ;
            bar_size = [fig_size(1), bar_width] ;
            color_map = reshape(color_map, [1, nb_color, 3]) ;
            X = list_pos ;
            Y = 0 ;
            
        case 'b'
            bar_pos = fig_pos ;
            bar_pos(2) = bar_pos(2)-bar_dist-bar_width ;
            bar_size = [fig_size(1), bar_width] ;
            color_map = reshape(color_map, [1, nb_color, 3]) ;
            X = list_pos ;
            Y = 0 ;
            
        case 'r'
            bar_pos = fig_pos ;
            bar_pos(1) = bar_pos(1)+fig_size(1)+bar_dist ;
            bar_size = [bar_width, fig_size(2)] ;
            color_map = reshape(color_map, [nb_color, 1, 3]) ;
            X = 0 ;
            Y = list_pos ;
            
        case 'l'
            bar_pos = fig_pos ;
            bar_pos(1) = bar_pos(1)-bar_dist-bar_width ;
            bar_size = [bar_width, fig_size(2)] ;
            color_map = reshape(color_map, [nb_color, 1, 3]) ;
            X = 0 ;
            Y = list_pos ;
        
        otherwise
            error([flag_pos, ' is an unknown flag...']) ;
    end
    axes_bar = axes('Units', 'pixels', 'Position', [bar_pos, bar_size]) ;
    
    %% Inserting the color bar
    imagesc(X, Y, color_map) ;
    set(axes_bar, 'YDir', 'normal') ;
    
    %% Inserting frame
    set(axes_bar,'visible','off') ;
    axes_bar.XTick = [] ;
    axes_bar.YTick = [] ;
    
    if LineWidth>0
        set(axes_bar,'visible','on') ;
        axes_bar.LineWidth = LineWidth ;
    end
    
    %% Inserting ticks
    if isempty(scale_func)
        scale_func = @(x)x ;
    end
    if ~isempty(struct_tick)
        % Number of ticks
        nb_tick = length(struct_tick.list) ;
        
        % Scaling between 0 and 1
        struct_tick.list = ...
            scale_law(struct_tick.list, val_min, val_max, scale_func) ;
        struct_tick.list = reshape(struct_tick.list, [nb_tick, 1]) ;
        
        % Turning the color into a list
        if size(struct_tick.color, 1) < nb_tick
            struct_tick.color = ...
                repmat(struct_tick.color(1,:), [nb_tick, 1]) ;
        end
        
        % Extracting position on the scale bars
        switch flag_pos
            case 't'
                xtick_pos = repmat(struct_tick.list, [1,2]) ;
                ytick_pos = [0.5*ones(nb_tick,1), ...
                    0.5*ones(nb_tick,1)-struct_tick.ratio] ;

            case 'b'
                xtick_pos = repmat(struct_tick.list, [1,2]) ;
                ytick_pos = [-0.5*ones(nb_tick,1), ...
                    -0.5*ones(nb_tick,1)+struct_tick.ratio] ;

            case 'r'
                xtick_pos = [+0.5*ones(nb_tick,1), ...
                    +0.5*ones(nb_tick,1)-struct_tick.ratio] ;
                ytick_pos = repmat(struct_tick.list, [1,2]) ;

            case 'l'
                xtick_pos = [-0.5*ones(nb_tick,1), ...
                    -0.5*ones(nb_tick,1)+struct_tick.ratio] ;
                ytick_pos = repmat(struct_tick.list, [1,2]) ;
        end
        
        % Inserting tick
        for tick = 1:nb_tick
            line(xtick_pos(tick,:), ytick_pos(tick,:), ...
                'LineWidth', struct_tick.LineWidth, ...
                'Color', struct_tick.color(tick,:)) ;
        end
        
    end
    
    %% Inserting tick labels
    if ~isempty(struct_tickLabel)
        % Number of ticks
        nb_tick = length(struct_tickLabel.list) ;
        
        % Scaling between 0 and 1
        struct_tickLabel.list = ...
            scale_law(struct_tickLabel.list, ...
            val_min, val_max, scale_func) ;
        struct_tickLabel.list = ...
            reshape(struct_tickLabel.list, [nb_tick, 1]) ;
        
        % Turning the angle into a list
        if ~isfield(struct_tickLabel, 'theta')
            struct_tickLabel.theta = 0 ;
        end
        if length(struct_tickLabel.theta) < nb_tick
            struct_tickLabel.theta = ...
                repmat(struct_tickLabel.theta(1), [nb_tick, 1]) ;
        end
        
        % Distance to the color bar
        if ~isfield(struct_tickLabel, 'dist')
            struct_tickLabel.dist = 0 ;
        end
        
        % Extracting position on the scale bars
        switch flag_pos
            case 't'
                xtick_pos = struct_tickLabel.list*fig_size(1) ;
                ytick_pos = zeros(nb_tick, 1) + ...
                    bar_width + struct_tickLabel.FontSize/2 + ...
                    struct_tickLabel.dist - 1;
                if ~isfield(struct_tickLabel, 'Alignment') % if empty
                    struct_tickLabel.Alignment = 'center' ;
                end

            case 'b'
                xtick_pos = struct_tickLabel.list*fig_size(1) ;
                ytick_pos = zeros(nb_tick, 1) - ...
                    struct_tickLabel.FontSize/2 - ...
                    struct_tickLabel.dist - 1  ;
                if ~isfield(struct_tickLabel, 'Alignment') % if empty
                    struct_tickLabel.Alignment = 'center' ;
                end

            case 'r'
                xtick_pos = zeros(nb_tick, 1) + ...
                    bar_width + ...
                    struct_tickLabel.dist - 1 ;
                ytick_pos = struct_tickLabel.list*fig_size(2) ;
                if ~isfield(struct_tickLabel, 'Alignment') % if empty
                    struct_tickLabel.Alignment = 'left' ;
                end

            case 'l'
                xtick_pos = zeros(nb_tick, 1) - ...
                    struct_tickLabel.dist ;
                ytick_pos = struct_tickLabel.list*fig_size(2)  - 1 ;
                if ~isfield(struct_tickLabel, 'Alignment') % if empty
                    struct_tickLabel.Alignment = 'right' ;
                end
        end
        
        
        % Turning the alignment into a list
        if ischar(struct_tickLabel.Alignment) % if not a cell
            struct_tickLabel.Alignment = {struct_tickLabel.Alignment} ;
        end
        if length(struct_tickLabel.Alignment) < nb_tick
            struct_tickLabel.Alignment = ...
                repmat({struct_tickLabel.Alignment{1}}, [nb_tick, 1]) ;
        end
        
        % Inserting tick
        for tick = 1:nb_tick
            h=text( ...
                xtick_pos(tick), ytick_pos(tick), ...
                struct_tickLabel.listLabel{tick}, ...
                'Units', 'pixels', ...
                'Interpreter', 'Latex');
            set(h, ...
                'Rotation', struct_tickLabel.theta(tick), ...
                'HorizontalAlignment', ...
                struct_tickLabel.Alignment{tick}, ...
                'FontSize', struct_tickLabel.FontSize);
        end
    end
    
    %% Going back to the main figure
    axes_main = axes('OuterPosition', [0, 0, 1, 1]) ;
    set(axes_main,'Units','pixels') ;
    set(axes_main, 'Position', pos);
    axis(axes_main, 'off') ;
    axis([lim_x, lim_y]) ;
    set(axes_main,'Units','normalized') ;
    
    pause(0.1) ;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to scale the law of a variable between 0 and 1
% 
% Created: 04/04/2019 (mm/dd/yyyy)
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #var# variable to scale
%
% #[val_min, val_max]# extremal value projected on 0 and 1
%
% #scale_func# scaling function between 0 and 1
%
%%%%%%%%
% Ouput
%%%%%%%%
% #scale_var# scaled variable between 0 and 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scale_var = scale_law(var, val_min, val_max, scale_func)
    
    %% Applying the scaling law
    scale_var = scale_func(var) ;
    val_max = scale_func(val_max) ;
    val_min = scale_func(val_min) ;

    %% Projecting between 0 and 1
    scale_var = (scale_var-val_min)/(val_max-val_min) ;
    scale_var = min(max(0,scale_var),1) ;
end