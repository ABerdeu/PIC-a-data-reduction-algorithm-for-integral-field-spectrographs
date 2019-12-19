%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get the exposure time of a fits file
% 
% Created: 01/26/2019 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #fits_name# name of the fits file
%
%%%%%%%%
% Ouput
%%%%%%%%
% #exposure# The exposure time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [exposure] = get_fits_exposure(fits_name)
    
    %% Fits info
    fits_info = fitsinfo(fits_name) ;
    
    %% Finding exposure time
    % ESO DET SEQ1 REALDIT
    exposure = str2double(get_fits_info_field(fits_info, ...
        'ESO DET SEQ1 DIT    =', ' ')) ;
end