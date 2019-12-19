%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to compute the cross-correlation of two pictures to find the
% relative positions
% 
% Created: 10/14/2015 (mm/dd/yyyy)
% Author: Anthony Berdeu (PhD student at CEA)
% Modified: 06/28/2018 (mm/dd/yyyy)
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #pic# the picture in which the ROI must be extracted
%
% #ROI# the pattern of the region of interest to extract
%
%%%%%%%%
% Ouput
%%%%%%%%
% #[x_cor, y_cor]# the position of the ROI in pic
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_cor, y_cor] = cross_correlation(pic, ROI)

    %% Initialisation
    [nb_y, nb_x] = size(pic) ;
    [nb_y_cor, nb_x_cor] = size(ROI) ;
    cor_pat = ROI-mean(ROI(:)) ;

    % Fourier transform of ROI
    fft_cor_pat = zeros(nb_y, nb_x) ;
    fft_cor_pat((1:nb_y_cor)+floor(nb_y/2)-floor(nb_y_cor/2), ...
        (1:nb_x_cor)+floor(nb_x/2)-floor(nb_x_cor/2)) = cor_pat ;
    fft_cor_pat = fftshift(fft2(ifftshift(fft_cor_pat))) ;

    % Fourier transform of pic
    fft_cor_pic = fftshift(fft2(ifftshift(pic-mean(pic(:))))) ;

    % Correlation in the Fourier domain
    cor = fft_cor_pat.*conj(fft_cor_pic) ;
    cor = fftshift(ifft2(ifftshift(cor))) ;
    cor = abs(cor) ;

    % Extracting center
    [y_cor_c, x_cor_c] = find(cor==max(cor(:))) ;
    x_cor = nb_x-x_cor_c ;
    y_cor = nb_y-y_cor_c ;
    
    % According to odd or even dimensions
    x_cor = x_cor+0.5+0.5*mod(nb_x_cor,2)+1-mod(nb_x,2) ;
    y_cor = y_cor+0.5+0.5*mod(nb_y_cor,2)+1-mod(nb_y,2) ;
end