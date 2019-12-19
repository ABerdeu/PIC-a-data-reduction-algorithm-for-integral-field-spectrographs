%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to convert a picture matrix to a RGB matrix for a given colormap
% and value range
% 
% Created: 04/24/2018 (mm/dd/yyyy)
% Modified: 09/18/2018 (mm/dd/yyyy) Minor corrections for the input
% default values
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
% Modified: 08/29/2019 (mm/dd/yyyy) Rename because insert_colorbar is now a
% function to add a colorbar in a Matlab figure
% Author: Anthony Berdeu (Centre de Recherche en Astrophysique de
% Lyon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #pic_in# the gray matrix to convert to a picture with a colorbar
%
% #dyn = [val_min, val_max]# dynamics range of the values of the colorbar.
% If empty, automatic selection of the minimal and maximal values of the
% picture.
%
% #options# scalebar options
%   -> nb_pix_bar: width of the colorbar (default: 50)
%   -> nb_pix_delta: width of separation between the picture and the
%       colorbar (default: 25)
%   -> cmap: colormap of the picture (default: gray scale)
%   -> nb_stick: number of value to add to the colorbar (default: 10)
%   -> str_format: string format (default: '%0.2f')
%   -> FontSize: font size (default: 18)
%   -> BGcolor: Background color (default white: [1,1,1])
%   -> nb_edge: width of the edges of the colorbar If empty, 0.
%   -> TextColor: Text color (default black: [0,0,0])
%
%%%%%%%%
% Ouput
%%%%%%%%
% #pic_out# the RGB picture
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pic_out] = insert_colorbar_pic(pic_in, dyn, options)
    
    %% Initialization
    [nb_y, ~] = size(pic_in) ;
    
    % dyn
    if nargin < 2 || isempty(dyn)
        dyn = [min(pic_in(:)), max(pic_in(:))] ;
    end
    
    % options
    if nargin < 3
        options = [] ;
    end
    
    % nb_pix_bar
    if isfield(options, 'nb_pix_bar')
        nb_pix_bar = options.nb_pix_bar ;
    else
        nb_pix_bar = 50 ;
    end
    
    % nb_pix_delta
    if isfield(options, 'nb_pix_delta')
        nb_pix_delta = options.nb_pix_delta ;
    else
        nb_pix_delta = 25 ;
    end
    
    % colormap
    if isfield(options, 'colormap')
        cmap = options.colormap ;
    else
        cmap = colormap(gray) ;
    end
    nb_map = length(cmap) ;
    
    % nb_stick
    if isfield(options, 'nb_stick')
        nb_stick = options.nb_stick ;
    else
        nb_stick = 10 ;
    end
    
    % str_format
    if isfield(options, 'str_format')
        str_format = options.str_format ;
    else
        str_format = '%0.2f' ;
    end
    
    % FontSize
    if isfield(options, 'FontSize')
        FontSize = options.FontSize ;
    else
        FontSize = 18 ;
    end
    
    % BGcolor
    if isfield(options, 'BGcolor')
        BGcolor = options.BGcolor ;
    else
        BGcolor = [1,1,1] ;
    end
    
    % TextColor
    if isfield(options, 'TextColor')
        TextColor = options.TextColor ;
    else
        TextColor = [0,0,0] ;
    end
    
    % nb_edge
    if isfield(options, 'nb_edge')
        nb_edge = options.nb_edge ;
    else
        nb_edge = 0 ;
    end
    
    %% Rescale gray picture
    pic_in = (pic_in-dyn(1))/(dyn(2)-dyn(1))*(nb_map-1)+1 ;
    pic_in = min(max(round(pic_in),1),nb_map) ;
    pic_out = ind2rgb(pic_in,cmap) ;
    
    %% Adding colorbar
    scale_bar = linspace(nb_map, 0, nb_y)' ;
    scale_bar = min(max(floor(scale_bar),0),nb_map-1) ;
    scale_bar = repmat(scale_bar+1, [1,nb_pix_bar]) ;
    scale_bar = ind2rgb(scale_bar,cmap) ;
    
    % Edges of the colorbar
    if nb_edge>0
        scale_bar(1:nb_edge,:,:) = 0 ;
        scale_bar(:,1:nb_edge,:) = 0 ;
        scale_bar(:,end-(0:(nb_edge-1)),:) = 0 ;
        scale_bar(end-(0:(nb_edge-1)),:,:) = 0 ;
    end
    
    background = cat(3, BGcolor(1)*ones(nb_y, nb_pix_delta), ...
        BGcolor(2)*ones(nb_y, nb_pix_delta), ...
        BGcolor(3)*ones(nb_y, nb_pix_delta)) ;
    pic_out = cat(2,pic_out, background, scale_bar) ;
    [nb_y, nb_x, ~] = size(pic_out) ;
    
    %% Adding values
    if nb_stick>1
        % Values
        val_stick = linspace(dyn(2), dyn(1), nb_stick) ;
        text_str = cell(nb_stick,1) ;
        nb_str = 0 ;
        for s = 1:nb_stick
            text_str{s} = sprintf(str_format, val_stick(s)) ;
            nb_str = max(nb_str, length(text_str{s})) ;
        end
        
        % Positions of the values
        position = linspace(1, nb_y, nb_stick)' ;
        position = [nb_x + ...
            zeros(nb_stick,1), ... 
            position-FontSize/2] ;
        
        % Extension of the pictures to put the texts
        background = cat(3, BGcolor(1)*ones(floor(FontSize/2)+1, nb_x), ...
            BGcolor(2)*ones(floor(FontSize/2)+1, nb_x), ...
            BGcolor(3)*ones(floor(FontSize/2)+1, nb_x)) ;
        pic_out = cat(1, background, pic_out, background) ; % y-axis
        
        [nb_y, ~, ~] = size(pic_out) ;
        background = cat(3, BGcolor(1)*ones(nb_y, nb_str*FontSize), ...
            BGcolor(2)*ones(nb_y, nb_str*FontSize), ...
            BGcolor(3)*ones(nb_y, nb_str*FontSize)) ;
        pic_out = cat(2, pic_out, background) ; % x-axis
        
        pic_out = insertText(pic_out, position ,text_str, ...
            'FontSize', FontSize, 'BoxOpacity', 0, 'TextColor', TextColor);
        
        % Remove extra lines and columns
        [nb_y, nb_x, ~] = size(pic_out) ;
        nb_x_remove = 0 ;
        test_background = ...
            sum(pic_out(:,nb_x-nb_x_remove,1)==BGcolor(1)) + ...
            sum(pic_out(:,nb_x-nb_x_remove,2)==BGcolor(2)) + ...
            sum(pic_out(:,nb_x-nb_x_remove,3)==BGcolor(3)) ;
        while test_background == 3*nb_y
            nb_x_remove = nb_x_remove+1 ;
            test_background = ...
                sum(pic_out(:,nb_x-nb_x_remove,1)==BGcolor(1)) + ...
                sum(pic_out(:,nb_x-nb_x_remove,2)==BGcolor(2)) + ...
                sum(pic_out(:,nb_x-nb_x_remove,3)==BGcolor(3)) ;
        end
        pic_out = pic_out(:,1:(nb_x-nb_x_remove),:) ;
        [nb_y, nb_x, ~] = size(pic_out) ;
        
        nb_y_remove = 0 ;
        test_background = ...
            sum(pic_out(nb_y-nb_y_remove,:,1)==BGcolor(1)) + ...
            sum(pic_out(nb_y-nb_y_remove,:,2)==BGcolor(2)) + ...
            sum(pic_out(nb_y-nb_y_remove,:,3)==BGcolor(3)) ;
        while test_background == 3*nb_x
            nb_y_remove = nb_y_remove+1 ;
            test_background = ...
                sum(pic_out(nb_y-nb_y_remove,:,1)==BGcolor(1)) + ...
                sum(pic_out(nb_y-nb_y_remove,:,2)==BGcolor(2)) + ...
                sum(pic_out(nb_y-nb_y_remove,:,3)==BGcolor(3)) ;
        end
        pic_out = pic_out(1:(nb_y-nb_y_remove),:,:) ;
        [nb_y, nb_x, ~] = size(pic_out) ;
        
        nb_y_remove = 0 ;
        test_background = ...
            sum(pic_out(1+nb_y_remove,:,1)==BGcolor(1)) + ...
            sum(pic_out(1+nb_y_remove,:,2)==BGcolor(2)) + ...
            sum(pic_out(1+nb_y_remove,:,3)==BGcolor(3)) ;
        while test_background == 3*nb_x
            nb_y_remove = nb_y_remove+1 ;
            test_background = ...
                sum(pic_out(1+nb_y_remove,:,1)==BGcolor(1)) + ...
                sum(pic_out(1+nb_y_remove,:,2)==BGcolor(2)) + ...
                sum(pic_out(1+nb_y_remove,:,3)==BGcolor(3)) ;
        end
        pic_out = pic_out((1+nb_y_remove):nb_y,:,:) ;
    end
end
