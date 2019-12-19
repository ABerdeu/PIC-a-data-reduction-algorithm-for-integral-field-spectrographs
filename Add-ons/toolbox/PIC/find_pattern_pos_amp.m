%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to find the position and the amplitude of a given pattern whose
% parameters are known
% 
% Created: 08/31/2018 (mm/dd/yyyy)
% Modified: 02/06/2018 (mm/dd/yyyy) Accounting for bad pixels and
% background
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #pic# picture in which the pattern is extracted
%
% #par_in = [y_in, x_in, amp_in]# input position and amplitude
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
%   RL2_method -> method for the robust penalization
%   noise_model: model of the noise to scale the residues
%   var_0: value of the variance of the data at null flux for the Poisson
%   noise model
%   eta: ratio to scale the model in the corresponding unit for the
%   Poisson noise
%   All the needed optimization options
%
% #pattern_model# flags to specify the pattern to fit
%   pattern_model.flag_norm -> Normalized?
%   pattern_model.flag_profile -> Gaussian / Hexagon / Moffat
%
% #pattern_car# cell containing the list of the pattern characteristic 
% (positions, amplitude, offset...) (if a field is empty, it becomes a
% variable to estimate)
% 
% #par_pat# parameters of the pattern
% 
% #offset# offset of the pattern
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
% #[y_pat, x_pat]# position of the gaussian pattern
%
% #amp_pat# amplitude of the gaussian pattern
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y_pat, x_pat, amp_pat] = find_pattern_pos_amp(pic, par_in, ...
    rad_ROI, pix, option_opti, pattern_model, par_pat, offset, ...
    flag_fit, I_BG, flag_BP)

    %% Initialization
    % Type of fitting
    if nargin < 9
        flag_fit = 'pos_amp' ;
    end
    
    % Extraction of the pattern of the picture
    [ind_y, ind_x] = get_list_index( ...
        [par_in(1,1), par_in(2,1)], [rad_ROI, rad_ROI], pix) ;
    if isnumeric(pic)
        % pic is a matrix
        if isfield(pix, 'f')
            % Extraction of a given list of frame
            pic = pic(ind_y, ind_x, pix.f) ;
        else
            pic = pic(ind_y, ind_x, 1:pix.nb_f) ;
        end
        
        % Background extraction
        if nargin>9 && ~isempty(I_BG)
            I_BG = I_BG(ind_y, ind_x) ;
        else
            I_BG = 0 ;
        end
        
        % Bad pixels flag
        if nargin>10 && ~isempty(flag_BP)
            flag_BP = flag_BP(ind_y, ind_x) ;
        else
            flag_BP = ones([length(ind_y), length(ind_x)]) ;
        end
        
    else
        % pic is a memmapfile
        if isfield(pix, 'f')
            % Extraction of a given list of frame
            pic = extract_bin(pic, {ind_y, ind_x, pix.f}, ...
                [pix.nb_y, pix.nb_x, pix.nb_f]) ;
        else
            pic = extract_bin(pic, {ind_y, ind_x, 1:pix.nb_f}, ...
                [pix.nb_y, pix.nb_x, pix.nb_f]) ;
        end
        
        % Background extraction
        if nargin>9 && ~isempty(I_BG)
            I_BG = extract_bin(I_BG, {ind_y, ind_x}, ...
                [pix.nb_y, pix.nb_x]) ;
        else
            I_BG = 0 ;
        end
        
        % Bad pixels flag
        if nargin>10 && ~isempty(flag_BP)
            flag_BP = extract_bin(flag_BP, {ind_y, ind_x}, ...
                [pix.nb_y, pix.nb_x]) ;
        else
            flag_BP = ones([length(ind_y), length(ind_x)]) ;
        end
    end
    
    % Extracted positions
    list_y = ind_y - floor(pix.nb_y/2+1) ;
    list_x = ind_x - floor(pix.nb_x/2+1) ;

    %% Taking into account the noise and the bad pixels in the scaling of
    % the residues
    flag_BP(flag_BP==0) = Inf ;
    if isfield(option_opti, 'noise_model') && ...
            strcmp(option_opti.noise_model, 'Poisson') % The Poisson noise
            % is dynamic
        option_opti.flag_s = flag_BP.*option_opti.flag_s ;
        option_opti.var_0 = option_opti.var_0+option_opti.eta*I_BG ;
        option_opti.noise_model = 'Poisson' ;
    else % The Poisson noise is in the scaling of the residues /!\ Only the
            % contribution of the background is taken into account...
        option_opti.flag_s = flag_BP.*option_opti.flag_s.* ...
            (option_opti.eta.*I_BG+option_opti.var_0).^0.5 ;
        option_opti.noise_model = 'none' ;
    end
    
    %% Iterative fitting
    switch flag_fit
        case 'pos_amp'
            % Constraints
            if ~strcmp('VMLMB', option_opti.method)
                const_min = [] ;
                const_max = [] ;
            else
                const_min = zeros(3, pix.nb_f) ;
                const_min(1) = par_in(1)-rad_ROI/2 ;
                const_min(2) = par_in(2)-rad_ROI/2 ;
                const_min(3,:) = 0 ;
                const_max = zeros(size(par_in)) ;
                const_max(1) = par_in(1)+rad_ROI/2 ;
                const_max(2) = par_in(2)+rad_ROI/2 ;
                const_max(3,:) = +Inf ;
                
            end
            
            % Building input for the pattern
            pattern = get_pattern([], [], [], offset, ...
                pattern_model.flag_profile, par_pat) ;
            
            % Cost
            option_opti_robust = option_opti ;
            option_opti_robust.method = option_opti_robust.RL2_method ;
            C = get_cost_pattern(pic, list_y, list_x, ...
                pattern_model, pattern, option_opti_robust, ...
                [true, true, false]) ;
            
            % Optimization
            par_it = par_in ;
            par_it = run_Opti(C, const_min, const_max, ...
                option_opti, par_it) ;
            par_in(1:3, 1:pix.nb_f) = par_it ;
        case 'pos'
            % Constraints
            if ~strcmp('VMLMB', option_opti.method)
                const_min = [] ;
                const_max = [] ;
            else
                const_min = zeros(2,1) ;
                const_max = zeros(2,1) ;
                const_min(1) = par_in(1)-rad_ROI/2 ;
                const_min(2) = par_in(2)-rad_ROI/2 ;
                const_max(1) = par_in(1)+rad_ROI/2 ;
                const_max(2) = par_in(2)+rad_ROI/2 ;
            end
            
            % Building input for the pattern
            pattern = get_pattern([], [], par_in(3, 1:pix.nb_f), ...
                offset, pattern_model.flag_profile, par_pat) ;
            
            % Cost
            option_opti_robust = option_opti ;
            option_opti_robust.method = option_opti_robust.RL2_method ;
            C = get_cost_pattern(pic, list_y, list_x, ...
                pattern_model, pattern, option_opti_robust, true) ;
            
            % Optimization
            par_it = par_in(1:2, 1) ;
            par_it = run_Opti(C, const_min, const_max, ...
                option_opti, par_it) ;
            par_in(1:2, 1) = par_it ;
        case 'amp'
            % Constraints
            if ~strcmp('VMLMB', option_opti.method)
                const_min = [] ;
                const_max = [] ;
            else
                const_min = zeros(1,pix.nb_f) ;
                const_max = +Inf*ones(1,pix.nb_f) ;
            end
            
            % Building input for the pattern
            pattern = get_pattern(par_in(1,1), par_in(2,1), [], ...
                offset, pattern_model.flag_profile, par_pat) ;
            
            % Cost
            option_opti_robust = option_opti ;
            option_opti_robust.method = option_opti_robust.RL2_method ;
            C = get_cost_pattern(pic, list_y, list_x, ...
                pattern_model, pattern, option_opti_robust, false) ;
            
            % Optimization
            par_it = par_in(3, 1:pix.nb_f) ;
            par_it = run_Opti(C, const_min, const_max, ...
                option_opti, par_it) ;
            par_in(3, 1:pix.nb_f) = par_it ;
        otherwise
            error(['Fitting option ''', flag_fit, ...
                ''' does not exist']) ;
    end

    y_pat = par_in(1, 1) ;
    x_pat = par_in(2, 1) ;
    amp_pat = par_in(3, :) ;
end