%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get the neighboring hexagons of a given point from the index
% map
% 
% Created: 05/18/2018 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #map_hex# map of the index of the hexagon in the picture
%   -> 0: the pixel is not associated with a known hexagon
%   -> i: index of the closest hexagon
%
% #pos = [y, x]# position of the point
%
% #side# side of the hexagon
%
% #theta# orientation of the hexagon
%
% #pix# pixel chracteristic of map_hex (default: dx=dy=1)
%   pix.nb_x -> number of the pixel on the x-axis
%   pix.nb_y -> number of the pixel on the y-axis
%   pix.dx -> size of the pixel on the x-axis
%   pix.dy -> size of the pixel on the y-axis
%
%%%%%%%%
% Ouput
%%%%%%%%
% #list_neigh# list of the neighboring hexagons
%   -> 0: the pixel is not associated with a known hexagon
%   -> i: index of the closest hexagon
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function list_neigh = get_neigh_hex(map_hex, pos, side, theta, pix)

    %% Initialization
    if nargin < 5
        [pix.nb_y, pix.nb_x] = size(map_hex) ;
        pix.dx = 1 ;
        pix.dy = 1 ;
    end

    %% Creating hexagon
    [x_hex, y_hex] = rot_2D(theta+(0:5)'*60, side, 0) ;
    x_hex = x_hex + pos(2) ;
    y_hex = y_hex + pos(1) ;
    ind_hex = pos2ind([y_hex, x_hex], pix) ;
    if size(ind_hex, 1) == 1
        list_neigh = map_hex(ind_hex(1)+pix.nb_x*(ind_hex(2)-1)) ;
    else
        list_neigh = map_hex(ind_hex(:,1)+pix.nb_x*(ind_hex(:,2)-1)) ;
    end
end