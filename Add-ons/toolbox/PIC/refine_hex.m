%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to refine the center of a hexagon knowing its neighbors
% 
% Created: 05/24/2018 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #pos_hex = [y_hex, x_hex]# initial center for the hexagon
%
% #list_c# the list of the known corners
%
% #list_pos = [y_list, x_list]# the list of the positions of the known
%   corners
%
% #side_hex# initial side of the hexagon
%
% #theta_hex# initial orientation of the hexagon
%
%%%%%%%%
% Ouput
%%%%%%%%
% #pos_out# the refined position of the hexagon
%
% #side_out# the refined side of the hexagon
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pos_out, side_out] = refine_hex(pos_hex, ...
    list_c, list_pos, side_hex, theta_hex)

    options = optimset('Display','none') ;
    delta_par = fminsearch( ...
        @(delta_par)get_cost_hex(pos_hex, list_c, list_pos, side_hex, ...
    theta_hex, delta_par), [0.1, 0.1, 0.05*side_hex], options) ;
    pos_out = [pos_hex(1) + delta_par(1), pos_hex(2) + delta_par(2)] ;
    side_out = side_hex + delta_par(3) ;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #pos_hex = [y_hex, x_hex]# initial center for the hexagon
%
% #list_c# the list of the known corners
%
% #list_pos = [y_list, x_list]# the list of the positions of the known
%   corners
%
% #side_hex# initial side of the hexagon
%
% #theta_hex# initial orientation of the hexagon
%
% #delta_par = [dy, dx, d_side]# current guess for the cost function
%
%%%%%%%%
% Ouput
%%%%%%%%
% #cost# cost of the update positions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cost = get_cost_hex(pos_hex, list_c, list_pos, side_hex, ...
    theta_hex, delta_par)

    %% Initialization
    d_y = delta_par(1) ;
    d_x = delta_par(2) ;
    d_side = delta_par(3) ;
    
    %% Corners positioning
    [x_c, y_c] = rot_2D(theta_hex+(list_c-1)*60, side_hex+d_side, 0) ;
    y_c = y_c+pos_hex(1)+d_y ;
    x_c = x_c+pos_hex(2)+d_x ;
    
    %% Errors computation
    cost = sqrt((list_pos(:,1)-y_c).^2+(list_pos(:,2)-x_c).^2) ;
    cost = sum(cost) / length(list_c) ;
    
% % % % % % % % % % % % % % % %     fig_2 = figure(2) ;
% % % % % % % % % % % % % % % %     pix_IFS.nb_x = 2048 ;
% % % % % % % % % % % % % % % %     pix_IFS.nb_y = 2048 ;
% % % % % % % % % % % % % % % %     pix_IFS.dy = 1 ;
% % % % % % % % % % % % % % % %     pix_IFS.dx = 1 ;
% % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % %     plot_rectangle(fig_2, ...
% % % % % % % % % % % % % % % %         pos_hex+[d_y, d_x], ...
% % % % % % % % % % % % % % % %         [5, 5], 0, ...
% % % % % % % % % % % % % % % %         pix_IFS, [1,1,0]) ;
% % % % % % % % % % % % % % % %     plot_rectangle(fig_2, ...
% % % % % % % % % % % % % % % %         [y_c, x_c], ...
% % % % % % % % % % % % % % % %         [5, 5], 0, ...
% % % % % % % % % % % % % % % %         pix_IFS, [1,0,0]) ;
% % % % % % % % % % % % % % % %     plot_circle(fig_2, ...
% % % % % % % % % % % % % % % %         list_pos, ...
% % % % % % % % % % % % % % % %         2.5, pix_IFS, [1,0,0]) ;
% % % % % % % % % % % % % % % %     pause(0.01) ;
    
end