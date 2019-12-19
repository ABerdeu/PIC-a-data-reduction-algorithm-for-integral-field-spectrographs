%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to insert a scale bar in a picture
% 
% Created: 08/05/2018 (mm/dd/yyyy)
% Author: Anthony Berdeu (Centre de Recherche en Astrophysique de
% Lyon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #h_fig# the figure handler
%
% #text_bar# text of the scale bar
%
% #width# width of the scale bar
%
% #d_pix# pixel pitch (in the scale bar unit
%
% #pos# the position of the center of scale bar
%
% #LineWidth# line width to of the scale bar
%
% #FontSize# font size of the text
%
% #struct_par# additionnal parameters
%   .Alignment  -> the text alignment (default: center)
%       {'left', 'center', 'right'}
%   .theta      -> the orientation (default: 0)
%   .dist       -> distance of the text to the scale bar (default: 0)
%   .color      -> color of the scale bar (default: [0,0,0])
%
%%%%%%%%
% Ouput
%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function insert_scalebar(h_fig, text_bar, width, d_pix, ...
    pos, LineWidth, FontSize, struct_par)
    
    %% Selecting figure
    figure(h_fig) ;
    axes_main = gca(h_fig) ;
    
    %% Initialization of the parameters
    if nargin < 8 || isempty(struct_par)
        struct_par = [] ;
    end
    if ~isfield(struct_par, 'Alignment')
        struct_par.Alignment = 'center' ;
    end
    if ~isfield(struct_par, 'theta')
        struct_par.theta = 0 ;
    end
    if ~isfield(struct_par, 'dist')
        struct_par.dist = 0 ;
    end
    if ~isfield(struct_par, 'color')
        struct_par.color = [0, 0, 0] ;
    end
    
    
    %% Inserting the scale bar
    list_pos = zeros(2,2) ;
    list_pos(1,1) = pos(1) - cosd(struct_par.theta)*width/(2*d_pix) ;
    list_pos(1,2) = pos(2) - sind(struct_par.theta)*width/(2*d_pix) ;
    list_pos(2,1) = pos(1) + cosd(struct_par.theta)*width/(2*d_pix) ;
    list_pos(2,2) = pos(2) + sind(struct_par.theta)*width/(2*d_pix) ;
    
    line(list_pos(:,1), list_pos(:,2), ...
        'LineWidth', LineWidth, ...
        'Color', struct_par.color) ;
    
    %% Inserting the text
    list_pos = zeros(1,2) ;
    list_pos(1) = pos(1) - sind(struct_par.theta)*struct_par.dist ;
    list_pos(2) = pos(2) + cosd(struct_par.theta)*struct_par.dist ;
    
    h = text(list_pos(1), list_pos(2), ...
        text_bar, 'Interpreter', 'Latex');
    set(h, ...
        'Rotation', struct_par.theta, ...
        'HorizontalAlignment', struct_par.Alignment, ...
        'FontSize', FontSize, ...
        'Color', struct_par.color);
    
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