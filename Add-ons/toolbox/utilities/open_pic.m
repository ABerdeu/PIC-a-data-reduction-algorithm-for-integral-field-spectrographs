%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to open a given picture by automatically determining the format
% 
% Created: 06/06/2018 (mm/dd/yyyy)
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #name_data# total name of the file
%
% #flip_y# flipping the y-axis for tif stack opening (default 0)
%   1 -> yes
%   0 -> no
%
% #flip_z# flipping the z-axis for tif stack opening (default 0)
%   1 -> yes
%   0 -> no
%
% #verbose_warning# display warning if the opened picture is not in double
%   format?
%
%%%%%%%%
% Ouput
%%%%%%%%
% #pic# the open picture
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pic = open_pic(name_data, flip_y, flip_z, verbose_warning)

    %% Initialization
    % Checking inputs
    if nargin < 2
        flip_y = 0 ;
    elseif isempty(flip_y)
        flip_y = 0 ;
    end
    if nargin < 3
        flip_z = 0 ;
    elseif isempty(flip_y)
        flip_z = 0 ;
    end
    if nargin < 4
        verbose_warning = false ;
    elseif isempty(verbose_warning)
        verbose_warning = false ;
    end
    
    % Get file extension
    [~ ,~ ,file_ext] = fileparts(name_data) ;
    
    %% Opening file
    switch file_ext
        case {'.png', '.jpg'}
            pic = imread(name_data) ;
        case '.fits'
            pic = fitsread(name_data) ;
        case {'.tif', '.tiff'}
            pic = open_tif([], name_data, flip_y, flip_z) ;
        otherwise
            error(['''', file_ext, ...
                ''' file extension is empty or unknown...']) ;
    end
    
    %% Converting to double if necessary
    who_pic = whos('pic') ;
    if ~strcmp('double', who_pic.class) 
        if verbose_warning
            warning('Converting opened picture to ''double''...') ;
        end
        pic = double(pic) ;
    end
end