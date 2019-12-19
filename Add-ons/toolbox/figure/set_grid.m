%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to set the grid of a figure
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
% #LineStyle# the line style of the grid
%
% #Color# the color
%
% #Alpha# the transparency
%
% #Layer# the layer
%   'top'
%   'bottom'
%
%%%%%%%%
% Ouput
%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function set_grid(h_fig, LineStyle, Color, Alpha, Layer)
    
    % Selecting figure
    figure(h_fig) ;

    % Getting axes
    axes = gca ;
    
    % Set the parameters
    axes.XGrid = 'on' ;
    axes.YGrid = 'on' ;
    axes.GridLineStyle = LineStyle ;
    axes.GridColor = Color ;
    axes.GridAlpha = Alpha ;
    axes.Layer = Layer ;
    
    

end