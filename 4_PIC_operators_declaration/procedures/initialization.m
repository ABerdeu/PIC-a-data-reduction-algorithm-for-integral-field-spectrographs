%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Procedure to initialize the simulation
% 
% Created: 05/31/2018 (mm/dd/yyyy)
% Created: 08/08/2018 (mm/dd/yyyy) New functions
% Modified: 10/25/2018 (mm/dd/yyyy) Taking into account the spectral flat
% fitting of the parmeters
% Created: 11/27/2018 (mm/dd/yyyy) Adaptation to polynomial laws
% Created: 03/13/2019 (mm/dd/yyyy) Adaptation to local polynomial laws and
% autocalibration on the data
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Initialization...') ;
tic_aux = tic ;

if isempty(pix.nb_x_pad)
    pix.nb_x_pad = pix.nb_x ;
end
if isempty(pix.nb_y_pad)
    pix.nb_y_pad = pix.nb_y ;
end

if c_cal>0
    suf_name = ['_it_', num2str(c_cal)] ;
else
    suf_name = '' ;
end

%% Opening data
% Coefficients of the polynomial fittings
coef_pol_dif = load([path_calib_spec, ...
    'coef_pol_dif', suf_name, '.txt'], '-ascii') ;
coef_pol_x = load([path_calib_spec, ...
    'coef_pol_x', suf_name, '.txt'], '-ascii') ;
coef_pol_y = load([path_calib_spec, ...
    'coef_pol_y', suf_name, '.txt'], '-ascii') ;
list_lambda = load([path_calib_spec, 'list_lambda.txt'], '-ascii') ;
lambda_0 = load([path_calib_spec, 'lambda_0.txt'], '-ascii') ;
trans = load([path_calib_spec, 'trans.txt'], '-ascii') ;
lamp_spec = load([path_calib_spec, 'lamp_spec.txt'], '-ascii') ;

% Number of lenslets
nb_lenslet = size(coef_pol_x, 1) ;

% Number of wavelength
nb_lambda = length(lamp_spec) ;

% Lenslet characteristic in the absolute grid
lenslet_side = ...
    (coef_pol_x(2:7,1) - coef_pol_x(1,1)) + ...
    (coef_pol_y(2:7,1) - coef_pol_y(1,1))*1i ;
lenslet_side = mean( abs(lenslet_side).* ...
    exp(1i*mod(angle(lenslet_side)*180/pi,60)*pi/180)) ;
lenslet_theta = 30 + ... 
    mod(180/pi*angle(lenslet_side)+30,60)-30 ;
lenslet_dist = abs(lenslet_side) ; % Distance between the lenslets
lenslet_side = lenslet_dist/sqrt(3) ;

% Fill factor
lenslet_side = lenslet_fill_factor/100 * ...
    lenslet_side ;

%% Spectra parameters
spectra = [] ;
spectra.trans = trans ;
spectra.list_lambda = list_lambda-lambda_0 ;
spectra.PSF_nb_x = 2*rad_ROI+1 ;
spectra.PSF_nb_y = 2*rad_ROI+1 ;
spectra.coef_pol_y = coef_pol_y ;
spectra.coef_pol_x = coef_pol_x ;
spectra.coef_pol_dif = coef_pol_dif ;


%% Selection of the lenslet in the field
list_x = coef_pol_x(: ,1) ;
list_y = coef_pol_y(: ,1) ;
