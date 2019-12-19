%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to initialize the excution
% 
% Created: 01/26/2019 (mm/dd/yyyy)
% Modified: 02/25/2019 (mm/dd/yyyy) Analyzing all dark current
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Opening flat data
disp('Opening flats...') ;
list_flat = dir([data.path, '*', data.flat_name, '*']) ;
nb_flat = length(list_flat) ;
IFS_flat_list = cell(nb_flat, 1) ;
exposure_list_flat = zeros(nb_flat, 2) ;
for flat = 1:nb_flat
    % Opening data
    IFS_flat_list{flat} = fitsread([data.path, list_flat(flat).name]) ;
    
    % Opening header
    fits_info = fitsinfo([data.path, list_flat(flat).name]) ;
    
    % ESO DET SEQ1 REALDIT
    exposure_list_flat(flat,1) = ...
        str2double(get_fits_info_field(fits_info, ...
        'ESO DET SEQ1 DIT    =', ' ')) ;
    
    LAMP = get_fits_info_field(fits_info, 'ESO INS2 CAL', '''') ;
    switch LAMP
        case 'OFF     '
            exposure_list_flat(flat,2) = 0 ;
            
        case 'BBL     '
            exposure_list_flat(flat,2) = 5 ;
            
        case 'HL1     '
            exposure_list_flat(flat,2) = 1 ;
            
        case 'HL2     '
            exposure_list_flat(flat,2) = 2 ;
            
        case 'HL3     '
            exposure_list_flat(flat,2) = 3 ;
            
        case 'HL4     '
            exposure_list_flat(flat,2) = 4 ;
            
        otherwise
            error([LAMP, 'in an unknown ESO INS2 CAL']) ;
    end
end

% Removing 'OFF' flat (done with lenslet array...)
IFS_flat_list = IFS_flat_list(exposure_list_flat(:,2)>0) ;
exposure_list_flat = exposure_list_flat(exposure_list_flat(:,2)>0,:) ;
nb_flat = length(IFS_flat_list) ;

disp('Opening flat: done!') ;



%% Opening dark background data
disp('Opening background darks...') ;
list_dark_BG = dir([data.path, '*', data.dark_BG_name, '*']) ;
nb_dark_BG = length(list_dark_BG) ;
IFS_dark_list_BG = cell(nb_dark_BG, 1) ;
exposure_list_dark_BG = zeros(nb_dark_BG,1) ;
for dark = 1:nb_dark_BG
    % Opening data
    IFS_dark_list_BG{dark} = fitsread([data.path, list_dark_BG(dark).name]) ;
    
    % Opening header
    fits_info = fitsinfo([data.path, list_dark_BG(dark).name]) ;
    
    % ESO DET SEQ1 REALDIT
    exposure_list_dark_BG(dark,1) = ...
        str2double(get_fits_info_field(fits_info, ...
        'ESO DET SEQ1 DIT    =', ' ')) ;
end
disp('Opening background darks: done!') ;


%% Parameter
pix_IFS.nb_y = size(IFS_flat_list{1},1) ;
pix_IFS.nb_x = size(IFS_flat_list{1},2) ;
pix_IFS.dy = 1 ;
pix_IFS.dx = 1 ;

