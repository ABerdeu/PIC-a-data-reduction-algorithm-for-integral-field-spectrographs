%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to declare a figure
% 
% Created: 07/24/2018 (mm/dd/yyyy)
% Author: Anthony Berdeu (Centre de Recherche en Astrophysique de
% Lyon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #fig_n# the figure number
%
% #fig_pos# the position of the figure on the screen (pixels)
%
% #fig_size# the size of the figure on the screen (pixels)
%
%%%%%%%%
% Ouput
%%%%%%%%
% #fig_h# the figure handler
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fig_h = declare_figure(fig_n, fig_pos, fig_size)
    
    fig_h = figure(fig_n) ;
    set(fig_h, 'rend', 'painters', 'pos', ...
        [fig_pos, fig_size]) ;
    clf(fig_h) ;
    pause(0.1) ;
end