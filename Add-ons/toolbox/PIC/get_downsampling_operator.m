%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get the parameters of a downsampling operation as well as its
% operator
% 
% Created: 09/05/2018 (mm/dd/yyyy)
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #rate_OS# oversampling rate (must be an integer)
%
% #[list_y_DS, list_x_DS]# List of the positons on which the downsampling
% is performed (assumed to have unitary increment) (if list_x_DS is not
% provided, the oversampling is mono-dimensional)
%
%%%%%%%%
% Ouput
%%%%%%%%
% #OpDS# Downsampling operator
%
% #[list_y_OS, list_x_OS]# Corresponding oversampling positions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OpDS, list_y_OS, list_x_OS] = get_downsampling_operator( ...
    rate_OS, list_y_DS, list_x_DS)

    if nargin<3 % 1D
        %% Oversampled positions
        list_y_OS = ((min(list_y_DS)-1):1/rate_OS:(max(list_y_DS)+1))' ;

        %% Building downsampling operation
        % Averaging value on the oversampled frame
        psf = zeros(length(list_y_OS), 1) ;
        psf(1:rate_OS, 1) = 1/rate_OS ;
        OpAvg = LinOpConv('PSF', psf, true, 1) ;

        % Down sampling by interpolation (offset to take into account that
        % the PSF shift the center...)
        Interp = LinOpInterp({list_y_OS}, ...
            list_y_DS+0.5*(rate_OS-1)/rate_OS) ;
    
        % Global down sampling operator
        OpDS = Interp*OpAvg ;
        
    else  % 2D
        %% Oversampled positions
        list_y_OS = ((min(list_y_DS)-1):1/rate_OS:(max(list_y_DS)+1))' ;
        list_x_OS = ((min(list_x_DS)-1):1/rate_OS:(max(list_x_DS)+1))' ;

        %% Building downsampling operation
        % Averaging value on the oversampled frame
        psf = zeros(length(list_y_OS), length(list_x_OS)) ;
        psf(1:rate_OS, 1:rate_OS) = 1/rate_OS.^2 ;
        OpAvg = LinOpConv('PSF', psf, true, 1:2) ;

        % List of downsample positions
        [xx, yy] = meshgrid(list_x_DS, list_y_DS) ;
        OpShape = LinOpShape(size(yy(:)), size(yy)) ;

        % Down sampling by interpolation (offset to take into account that
        % the PSF shift the center...)
        Interp = LinOpInterp({list_y_OS, list_x_OS}, ...
            [yy(:), xx(:)]+0.5*(rate_OS-1)/rate_OS) ;
    
        % Global down sampling operator
        OpDS = OpShape*Interp*OpAvg ;
    end
end