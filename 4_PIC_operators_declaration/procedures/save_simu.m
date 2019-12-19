%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure save the simulation
% 
% Created: 04/27/2018 (mm/dd/yy)
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Saving ...') ;

% Variables
save([save_path, 'lenslet_grid'], 'lenslet_grid') ;
save([save_path, 'pix'], 'pix') ;
save([save_path, 'sensor'], 'sensor') ;
save([save_path, 'plot_ind'], 'plot_ind') ;

% Central slice and the hexagonal grid
save_fits(hypercube(:,:,plot_ind), 'slice', save_path) ;

% Convolution kernel
save_fits(Conv_hexa_2D*test, 'Conv_hexa_2D', save_path) ;

% Convolution
save_fits(hypercube_conv(:,:,plot_ind), 'slice_conv', save_path) ;

% Interpolation
save_fits(pic_interp, 'slice_conv_interp', save_path) ;

% Projection
save_fits(data_sim, 'IFS_sim', save_path) ;


disp('Saving done!') ;