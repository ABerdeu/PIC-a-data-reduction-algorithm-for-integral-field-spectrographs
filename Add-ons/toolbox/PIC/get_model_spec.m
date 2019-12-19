%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get the model of the projection of a set of spectra. The
% elementary pattern is assumed to be gaussian and the offset is assumed to
% be null.
% 
% Created: 03/07/2019 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #[list_y, list_x]# cartesian positions on which the model is computed
%
% #[deg_pol_y, deg_pol_x, deg_pol_dif]# degree of the polynomial laws of
% the model for the yx-position and the diffraction model
% 
% #list_lambda# list of the wavelengths in the model
% 
% #nb_spec# number of spectra
%
% #flag_fit_spec# 
%   ->  'SPEC_SHARED'	-> fitting the spectrum (assumed to be identical to
%       all projected spectra) knowing the transmission of each spectrum
%   ->  'SPEC'          -> fitting a spectrum (for each projected spectrum)
%       knowing the transmission of each spectrum
%   ->  'TRANS'         -> fitting the transmission of each spectrum with a
%       knwon global spectrum (assumed to be identical to all projected
%       spectra)
%
% #list_spec#
%   ->  flag_fit_spec = true    -> list of the transmission of each 
%           spectrum
%   ->  flag_fit_spec = false   -> common projected spectrum
%
% #pattern_model# flags to specify the pattern to fit
%   ->  pattern_model.flag_norm -> Normalized?
%   ->  pattern_model.flag_profile -> Gaussian / Hexagon
%   ->  pattern_model.oversampling -> oversampling of the pattern 
%   	(default: 1)
%%%%%%%%
% Ouput
%%%%%%%%
% #Spec_model# Global forward model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Spec_model] = get_model_spec(list_y, list_x, ...
    deg_pol_y, deg_pol_x, deg_pol_dif, list_lambda, nb_spec, ...
    flag_fit_spec, list_spec, pattern_model)


    %% Initialization
    % Number of coefficients in the polynomial laws
    nb_coef_y = deg_pol_y+1 ;
    nb_coef_x = deg_pol_x+1 ;
    nb_coef_dif = deg_pol_dif+1 ;
    
    % Number of wavelength
    nb_lambda = length(list_lambda) ;
    switch flag_fit_spec
        case {'SPEC', 'SPEC_SHARED'}
            % Fitting a spectrum for each spectrum
            coef_amp = reshape(list_spec, [nb_spec, 1]) ;
            coef_amp = repmat(coef_amp, [1, nb_lambda]) ;
            nb_amp = nb_lambda ;
            
        case 'TRANS' % Fitting the individual transmission
            coef_amp = reshape(list_spec, [1, nb_lambda]) ;
            coef_amp = repmat(coef_amp, [nb_spec, 1]) ;
            nb_amp = nb_spec ;
            
        otherwise
            error([flag_fit_spec, ' is a unknown flag.']) ;
        
    end
    
    % y-law polynomial
    Pol_y = zeros(nb_lambda, nb_coef_y) ;
    for d = 0:deg_pol_y
        Pol_y(:,d+1) = list_lambda.^d ;
    end

    % x-law polynomial
    Pol_x = zeros(nb_lambda, nb_coef_x) ;
    for d = 0:deg_pol_x
        Pol_x(:,d+1) = list_lambda.^d ;
    end

    % Diffraction polynomial
    Pol_dif = zeros(nb_lambda, nb_coef_dif) ;
    for d = 0:deg_pol_dif
        Pol_dif(:,d+1) = list_lambda.^d ;
    end
    
    %% Loop on the spectra
    list_mod_spec = cell(nb_spec, 1) ;
    for spec = 1:nb_spec
        % Loop on the wavelength
        list_mod_lambda = cell(nb_lambda, 1) ;
        for lambda = 1:nb_lambda
            % Matrix to select the coefficient of the spectrum
            switch flag_fit_spec
                case 'SPEC'
                    M_sel = zeros(4, ...
                        (nb_coef_y + nb_coef_x + nb_coef_dif + ...
                        nb_amp) * nb_spec) ;
                    
                case 'TRANS'
                    M_sel = zeros(4, ...
                        (nb_coef_y + nb_coef_x + nb_coef_dif) * ...
                        nb_spec + nb_spec) ;

                case 'SPEC_SHARED'
                    M_sel = zeros(4, ...
                        (nb_coef_y + nb_coef_x + nb_coef_dif) * ...
                        nb_spec + nb_amp) ;
            end
            
            % Selection of the coefficient of the y-law for the spec-th
            % spectrum
            M_sel(1, (spec-1)*nb_coef_y + (1:nb_coef_y)) = ...
                Pol_y(lambda,:) ;
            
            % Selection of the coefficient of the x-law for the spec-th
            % spectrum
            M_sel(2, nb_coef_y*nb_spec + ...
                (spec-1)*nb_coef_x + (1:nb_coef_x)) = Pol_x(lambda,:) ;
            
            % Selection of the coefficient of the dif-law for the spec-th
            % spectrum
            M_sel(3, (nb_coef_y+nb_coef_x)*nb_spec + ...
                (spec-1)*nb_coef_dif + (1:nb_coef_dif)) = ...
                Pol_dif(lambda,:) ;
            
            % Selection of the amplitude
            switch flag_fit_spec
                case 'SPEC'
                    M_sel(4, ...
                        (nb_coef_y+nb_coef_x+nb_coef_dif)*nb_spec + ...
                        nb_amp*(spec-1)+lambda) = coef_amp(spec, lambda) ;
                    
                case'TRANS'
                    M_sel(4, ...
                        (nb_coef_y+nb_coef_x+nb_coef_dif)*nb_spec + ...
                        +spec) = coef_amp(spec, lambda) ;

                case 'SPEC_SHARED'
                    M_sel(4, ...
                        (nb_coef_y+nb_coef_x+nb_coef_dif)*nb_spec + ...
                        lambda) = coef_amp(spec, lambda) ;
            end
            
            % Computation of the elementary PSF
            list_mod_lambda{lambda} = get_pattern_operator( ...
                list_y, list_x, [], [], [], [], 0, pattern_model) * ...
                LinOpMatrix(M_sel) ;
        end
        
        % Forward model for the spectrum
        list_mod_spec{spec} = MapSummation(list_mod_lambda, ...
            ones(nb_lambda,1)) ;
    end
    
    %% Global forward model
    Spec_model = MapSummation(list_mod_spec, ones(nb_spec,1)) ;
end