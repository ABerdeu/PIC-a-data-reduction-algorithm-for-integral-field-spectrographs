%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to find the parameters of a given pattern whose position and the
% amplitude parameters are known
% 
% Created: 08/31/2018 (mm/dd/yyyy)
% Modified: 02/06/2018 (mm/dd/yyyy) Accounting for bad pixels and
% background
% Modfied: 02/26/2019 (mm/dd/yyyy) Hexagon pattern
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #pic# picture in which the pattern is extracted
%
% #par_in# input parameters
% 
% #rad_ROI# radius of the region of interest to extract
%
% #pix# pixel chracteristic
%   pix.nb_x -> number of the pixel on the x-axis
%   pix.nb_y -> number of the pixel on the y-axis
%   pix.nb_f -> number of frames in the picture
%   pix.dx -> size of the pixel on the x-axis
%   pix.dy -> size of the pixel on the y-axis
%
% #option_opti# structure with the options to fit the pattern by an
% iterative matching
%   flag_s -> method to scale the residue before applying the
%       objective function (default: none)
%   RL2_method -> method for the reweighted least square
%   All the needed optimization options
%
% #pattern_model# flags to specify the pattern to fit
%   pattern_model.flag_norm -> Normalized?
%   pattern_model.flag_profile -> Gaussian / Hexagon / Moffat
%   pattern_model.scale -> the scale to apply at each wavelength (if
%   needed)
%
% #list_y, list_x# list of the positions of the the pattern's realization
% in the picture
% 
% #list_amp# of the amplitude of the pattern's realization in the picture
% at each frame (or is considered to be the same at all frames)
%
% #offset# offset of the patterns
%
% #flag_fit# flag to specify what is fitted
%   pos_amp -> fitting both position and amplitude (default)
%   pos -> fitting the position
%   amp -> fitting the amplitude
%
% #I_BG# Background intensity
%
% #flag_BP# flag on the bad pixels
%
%%%%%%%%
% Ouput
%%%%%%%%
% #par_out# output parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function par_out = find_pattern_par(pic, par_in, rad_ROI, pix, ...
    option_opti, pattern_model, list_y, list_x, list_amp, offset, ...
    I_BG, flag_BP)

    %% Initialization of the optimization
    if nargin<11
        I_BG = [] ;
    end
    if nargin<12
        flag_BP = [] ;
    end
    
    % Constraints
    switch pattern_model.flag_profile
        case {'Gaussian', 'Hexagon'}
            % Constraints
            if ~strcmp('VMLMB', option_opti.method)
                const_min = [] ;
                const_max = [] ;
            else
                const_min = 0.1 ;
                const_max = Inf ;
            end
        case 'Moffat'
            % Constraints
            if ~strcmp('VMLMB', option_opti.method)
                const_min = [] ;
                const_max = [] ;
            else
                const_min = [0.1; 1] ;
                const_max = [+Inf; +Inf] ;
            end
        otherwise
            error([pattern_model.flag_profile, ...
                ': Unkown model for the pattern...']) ;
    end
    
    % Number of parameters and wavelengths to fit
    [nb_par, nb_l] = size(par_in) ;
    nb_l_pic = size(pic, 4) ; % Taking in account that pic can be a map of 
    % residue at each wavelength
    
    % Fitting the parameters
    if nb_l > 1 && pattern_model.flag_scale % They are several wavelength
        % to commonly fit using a scaling factor
        for l = 1:nb_l
            % Extracting pictures for wavelength l
            pic_l = pic(:,:,:,min(l,nb_l_pic)) ;
            % Getting the cost for the wavelength
            C_l = get_cost_l(pic_l, rad_ROI, pix, option_opti, ...
                pattern_model, list_y(:,l), list_x(:,l), list_amp(:,l), ...
                offset, I_BG, flag_BP) ;
            
            % Scaling the scale parameter
            diag = ones(nb_par,1) ;
            diag(1) = pattern_model.scale(l) ;
            C_l = C_l*LinOpDiag([nb_par,1], diag) ;
            
            if l==1
                C = C_l ;
            else
                C = C + C_l ;
            end
        end
        
        % Global optimization
        par_out = run_Opti(C, const_min, const_max, ...
            option_opti, par_in(:,1)) ;
        
        % Scaling the output at each wavelength
        par_out = repmat(par_out, [1, nb_l]) ;
        for l = 1:nb_l
            par_out(1,l) = par_out(1,l)*pattern_model.scale(l) ;
        end
    else
        par_out = par_in ;
        for l = 1:nb_l
            % Extracting pictures for wavelength l
            pic_l = pic(:,:,:,min(l,nb_l_pic)) ;
            % Getting the cost for the wavelength
            C = get_cost_l(pic_l, rad_ROI, pix, option_opti, ...
                pattern_model, list_y(:,l), list_x(:,l), list_amp(:,l), ...
                offset, I_BG, flag_BP) ;
            
            % Optimization
            par_out(:,l) = run_Opti(C, const_min, const_max, ...
                option_opti, par_in(:,l)) ;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get the cost for a wavelength
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #pic# picture in which the pattern is extracted
%
% #rad_ROI# radius of the region of interest to extract
%
% #pix# pixel chracteristic
%   pix.nb_x -> number of the pixel on the x-axis
%   pix.nb_y -> number of the pixel on the y-axis
%   pix.nb_f -> number of frames in the picture
%   pix.dx -> size of the pixel on the x-axis
%   pix.dy -> size of the pixel on the y-axis
%
% #option_opti# structure with the options to fit the pattern by an
% iterative matching
%   flag_s -> method to scale the residue before applying the
%       objective function (default: none)
%   RL2_method -> method for the reweighted least square
%   All the needed optimization options
%
% #pattern_model# flags to specify the pattern to fit
%   pattern_model.flag_norm -> Normalized?
%   pattern_model.flag_profile -> Gaussian / Moffat
%
% #list_y, list_x# list of the positions of the the pattern's realization
% in the picture
% 
% #list_amp# of the amplitude of the pattern's realization in the picture
% at each frame (or is considered to be the same at all frames)
%
% #offset# offset of the patterns
%
% #flag_fit# flag to specify what is fitted
%   pos_amp -> fitting both position and amplitude (default)
%   pos -> fitting the position
%   amp -> fitting the amplitude
%
% #I_BG# Background intensity
%
% #flag_BP# flag on the bad pixels
%
%%%%%%%%
% Ouput
%%%%%%%%
% #C_l# the cost for the given wavelength
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C_l = get_cost_l(pic, rad_ROI, pix, option_opti, ...
    pattern_model, list_y, list_x, list_amp, offset, I_BG, flag_BP)
    
    % Number of pattern    
    nb_pat = length(list_y) ;
    for pat = 1:nb_pat
        % Extraction of the pattern of the picture
        [ind_y, ind_x] = get_list_index( ...
            [list_y(pat), list_x(pat)], [rad_ROI, rad_ROI], ...
            pix) ;
        if isnumeric(pic)
            % pic is a matrix
            if isfield(pix, 'f')
                % Extraction of a given list of frame
                pic_extract = pic(ind_y, ind_x, pix.f) ;
            else
                pic_extract = pic(ind_y, ind_x, 1:pix.nb_f) ;
            end

            % Background extraction
            if ~isempty(I_BG)
                I_BG_extract = I_BG(ind_y, ind_x) ;
            else
                I_BG_extract = 0 ;
            end

            % Bad pixels flag
            if ~isempty(flag_BP)
                flag_BP_extract = flag_BP(ind_y, ind_x) ;
            else
                flag_BP_extract = ones([length(ind_y), length(ind_x)]) ;
            end

        else
            % pic is a memmapfile
            if isfield(pix, 'f')
                % Extraction of a given list of frame
                pic_extract = extract_bin(pic, {ind_y, ind_x, pix.f}, ...
                    [pix.nb_y, pix.nb_x, pix.nb_f]) ;
            else
                pic_extract = extract_bin(pic, {ind_y, ind_x, 1:pix.nb_f}, ...
                    [pix.nb_y, pix.nb_x, pix.nb_f]) ;
            end

            % Background extraction
            if nargin>9 && ~isempty(I_BG)
                I_BG_extract = extract_bin(I_BG, {ind_y, ind_x}, ...
                    [pix.nb_y, pix.nb_x]) ;
            else
                I_BG_extract = 0 ;
            end

            % Bad pixels flag
            if nargin>10 && ~isempty(flag_BP)
                flag_BP_extract = extract_bin(flag_BP, {ind_y, ind_x}, ...
                    [pix.nb_y, pix.nb_x]) ;
            else
                flag_BP_extract = ones([length(ind_y), length(ind_x)]) ;
            end
        end
        
        % Taking into account the noise and the bad pixels in the scaling
        % of the residues
        flag_BP_extract(flag_BP_extract==0) = Inf ;
        if isfield(option_opti, 'noise_model') && ...
                strcmp(option_opti.noise_model, 'Poisson') % The Poisson
                % noise is dynamic
            option_opti.flag_s = flag_BP_extract.*option_opti.flag_s ;
            option_opti.var_0 = option_opti.var_0 + ...
                option_opti.eta*I_BG_extract ;
            option_opti.noise_model = 'Poisson' ;
        else % The Poisson noise is in the scaling of the residues /!\ Only
                % the contribution of the background is taken into
                % account...
            option_opti.flag_s = flag_BP_extract.*option_opti.flag_s.* ...
                (option_opti.eta.*I_BG_extract+option_opti.var_0).^0.5 ;
            option_opti.noise_model = 'none' ;
        end
        option_opti.method = option_opti.RL2_method ;
        
        % Building input for the pattern
        pattern = get_pattern(list_y(pat), list_x(pat), ...
            list_amp(pat,:), offset, pattern_model.flag_profile, []) ;

        % Cost
        pos_y = ind_y - floor(pix.nb_y/2+1) ;
        pos_x = ind_x - floor(pix.nb_x/2+1) ;
        C_aux = get_cost_pattern(pic_extract, pos_y, pos_x, ...
            pattern_model, pattern, option_opti) ;
        
        % Normalization
        C_aux = 1/numel(pic_extract)*C_aux ;
        
        if pat==1
            C_l = C_aux ;
        else
            C_l = C_l + C_aux ;
        end
    end
end