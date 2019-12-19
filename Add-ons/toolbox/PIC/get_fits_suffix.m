%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get the suffix to apply on a fits file
% 
% Created: 01/24/2019 (mm/dd/yyyy)
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
% #suffix# The wanted suffix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [suffix] = get_fits_suffix(fits_name)
    
    
    %% Fits info
    fits_info = fitsinfo(fits_name) ;
    
    %% Finding acquisition type
    % ESO DPR CATG
    CATG = get_fits_info_field(fits_info, ...
        'ESO DPR CATG', '''') ;
    
    % ESO DPR TYPE
    TYPE = get_fits_info_field(fits_info, ...
        'ESO DPR TYPE', '''') ;
    
    % ESO DPR TYPE
    PRIS = get_fits_info_field(fits_info, ...
        'ESO INS2 OPTI2 ID', '''') ;
    
    %% Building suffix
    switch CATG
        case 'CALIB   '
            %% Calibration file
            suffix = '_CALIB_' ;
            switch TYPE
                case 'WAVE,LAMP'
                    suffix = [suffix, 'WAVE_LAMP'] ;
                    
                case 'FLAT,LAMP'
                        LAMP = get_fits_info_field(fits_info, 'ESO INS2 CAL', '''') ;
                        switch LAMP
                            case 'OFF     '
                                suffix = [suffix, 'SPEC_FLAT_LAMP'] ;

                            case 'BBL     '
                                suffix = [suffix, 'FLAT_LAMP_BBL'] ;

                            case 'HL1     '
                                suffix = [suffix, 'FLAT_LAMP_HL1'] ;

                            case 'HL2     '
                                suffix = [suffix, 'FLAT_LAMP_HL2'] ;

                            case 'HL3     '
                                suffix = [suffix, 'FLAT_LAMP_HL3'] ;

                            case 'HL4     '
                                suffix = [suffix, 'FLAT_LAMP_HL4'] ;

                            otherwise
                                error([LAMP, ...
                                    'in an unknown ESO INS2 CAL']) ;
                        end
                    
                case 'DARK    '
                    suffix = [suffix, 'DARK'] ;
                    
                case 'DARK,BACKGROUND'
                    suffix = [suffix, 'DARK_BACKGROUND'] ;
                    
                case 'SPECPOS,LAMP'
                    suffix = [suffix, 'LAMP_SPECPOS'] ;
                    
                case 'LAMP,DISTORT'
                    suffix = [suffix, 'LAMP_DISTORT'] ;
                    
                case 'FLAT,LAMP,RONGAIN'
                    suffix = [suffix, 'FLAT_LAMP_RONGAIN'] ;
                    
                case 'STD,FLUX'
                    suffix = [suffix, 'FLUX_STD'] ;
                    
                otherwise
                    % Unknown file
                    error([TYPE, ' is an unknown ESO DPR TYPE']) ;
            end
                    
            
        case 'SCIENCE '
            %% Science file
            suffix = '_SCIENCE_' ;
            switch TYPE
                case 'OBJECT,CENTER'
                    suffix = [suffix, 'OBJECT_CENTER'] ;
                    
                case 'OBJECT  '
                    suffix = [suffix, 'OBJECT '] ;
                    
                case 'SKY     '
                    suffix = [suffix, 'SKY'] ;
                    
                case 'OBJECT,FLUX'
                    suffix = [suffix, 'OBJECT_CENTER_FLUX'] ;
                    
                otherwise
                    % Unknown file
                    error([TYPE, 'is an unknown ESO DPR TYPE']) ;
            end
            
        otherwise
            %% Unknown file
            error([CATG, 'is an unknown ESO DPR CATG']) ;
    end
    
    % Adding the prism
    switch PRIS
        case 'OPEN    '
            suffix = [suffix, '_OPEN'] ;

        otherwise
            suffix = [suffix, '_', PRIS] ;
    end
    
end