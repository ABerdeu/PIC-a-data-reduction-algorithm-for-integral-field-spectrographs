%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to analyze the spectral flat
% 
% Created: 02/25/2019 (mm/dd/yyyy)
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Opening spectral flat
IFS_spec_flat = dir([data.path, '*', data.spectral_flat_name, '*']) ;
if length(IFS_spec_flat) > 1
    warning(['There are several spectral flat in the specified', ...
        ' directory... Only the first one is taken into account...']) ;
end
IFS_spec_flat = [data.path, IFS_spec_flat(1).name] ;
IFS_spec_flat = fitsread(IFS_spec_flat) ;
IFS_spec_flat = mean(IFS_spec_flat,3) ; % Stacking the different frames

%% Parameter
pix_IFS.nb_y = size(IFS_spec_flat,1) ;
pix_IFS.nb_x = size(IFS_spec_flat,2) ;
pix_IFS.dy = 1 ;
pix_IFS.dx = 1 ;

%% Analyzing spectral flat
% Median filter (normalization by the median)
IFS_mask = IFS_spec_flat ;
IFS_mask(isnan(IFS_mask)) = 0 ;
IFS_mask = medfilt2(get_gray(IFS_mask, 2), [rad_med, rad_med]) ;

% Thresholding
IFS_mask = IFS_mask>threshold_data ;

% Closing picture
IFS_mask = imdilate(IFS_mask, strel('disk', rad_dil)) ;
% IFS_BG = imopen(IFS_BG, strel('disk', 2*rad_dil)) ;
% IFS_BG = imdilate(IFS_BG, strel('disk', rad_dil)) ;

% Removing edges
IFS_mask(1:nb_BP_edges, :) = false ;
IFS_mask(:, 1:nb_BP_edges) = false ;
IFS_mask(pix_IFS.nb_y+1+(-nb_BP_edges:-1), :) = false ;
IFS_mask(:, pix_IFS.nb_x+1+(-nb_BP_edges:-1)) = false ;

% Correcting acquisition
save_fits(uint8(IFS_mask*1), 'IFS_mask', save_path) ;
% save_fits(IFS_spec_flat, 'IFS_spec_flat', save_path) ;

%% Clearing
clear('IFS_spec_flat') ;

