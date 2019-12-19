%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to find the position and the amplitude of a gaussian pattern
% It performs a rough determination by a local maximization and then refine
% the position by fitting a gaussian pattern.
% 
% Created: 05/02/2018 (mm/dd/yyyy)
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
% #option_max# structure with the options of the localization with the
% maximal value (if empty, this step is not performed)
%   flag_method -> method to determine the local maximum
%   func -> function to apply to the picture (default: identity)
%   flag_conv -> several estimations until the position reaches convergence
%       (default: false)
%
% #option_opti# structure with the options to fit the gaussian position
% by an iterative matching (if empty, this step is not performed)
%   flag_s -> method to scale the residue before applying the
%       objective function (default: none)
%   RL2_method -> method for the reweighted least square
%   flag_conv -> several estimations until the position reaches convergence
%       (default: false)
%   VMLMB options
%
% #theta# orientation (in degree) of the gaussian pattern
% 
% #sig_par# elongation of the gaussian pattern along the
% direction given by theta
% 
% #sig_perp# elongation of the gaussian pattern along the
% perpendicular direction given by theta
% 
% #offset# offset of the gaussian pattern
%
% #flag_fit# flag to specify what is fitted
%   pos_amp -> fitting both position and amplitude (default)
%   pos -> fitting the position
%   amp -> fitting the amplitude
%
%%%%%%%%
% Ouput
%%%%%%%%
% #[y_gauss, x_gauss]# position of the gaussian pattern
%
% #amp_gauss# amplitude of the gaussian pattern
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y_gauss, x_gauss, amp_gauss] = find_gauss_pos_amp(pic, ...
    par_in, rad_ROI, pix, option_max, option_opti, ...
    theta, sig_par, sig_perp, offset, flag_fit)

    %% Initialization
    if nargin<11
        flag_fit = 'pos_amp' ;
    end
    
    %% Local maximum on the average picture
    if ~isempty(option_max)
        if ~strcmp(flag_fit, 'pos_amp')
            error(['Fitting option ''', flag_fit, ''' not possible ', ...
                ' with a predection by maximum fitting...', ...
                ' Please use ''pos_amp''']) ;
        end
        if ~isnumeric(pic)
            error(['Fitting option ''', flag_fit, ''' not possible ', ...
                ' with a non numeric picture...']) ;
        end
        [m_y, m_x, m] = find_local_maximum(mean(pic,3), ...
            [par_in(1,1), par_in(2,1)], rad_ROI, ...
            option_max.max_method, option_max.func, option_max.max_conv) ;
        par_in(1,1) = m_y ;
        par_in(2,1) = m_x ;
        par_in(3, 1:pix.nb_f) = m - offset ;
    end

    %% Iterative fitting
    if ~isempty(option_opti)
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
        end
        list_y = ind_y - floor(pix.nb_y/2+1) ;
        list_x = ind_x - floor(pix.nb_x/2+1) ;
        
        % Refining with a gaussian fitting
        switch flag_fit
            case 'pos_amp'
                % Constraints
                if ~strcmp('VMLMB', option_opti.method)
                    const_min = [] ;
                    const_max = [] ;
                else
                    const_min = zeros(3, pix.nb_f) ;
                    const_min(1,:) = -Inf ;
                    const_min(2,:) = -Inf ;
                    const_min(3,:) = 0 ;
                    const_max = zeros(size(par_in)) ;
                    const_max(1,:) = +Inf ;
                    const_max(2,:) = +Inf ;
                    const_max(3,:) = +Inf ;
                end
                % Cost
                C = get_cost_gauss(pic, list_y, list_x, ...
                    option_opti.flag_s, option_opti.RL2_method, ...
                    [], [], theta, sig_par, sig_perp, [], offset, ...
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
                    const_min(1) = -Inf ;
                    const_min(2) = -Inf ;
                    const_max = -const_min ;
                end
                C = get_cost_gauss(pic, list_y, list_x, ...
                    option_opti.flag_s, option_opti.RL2_method, ...
                    [], [], theta, sig_par, sig_perp, ...
                    par_in(3, 1:pix.nb_f), offset, true) ;
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
                % Cost
                C = get_cost_gauss(pic, list_y, list_x, ...
                    option_opti.flag_s, option_opti.RL2_method, ...
                    par_in(1,1), par_in(2,1), theta, sig_par, sig_perp, ...
                    [], offset, false) ;
                % Optimization
                par_it = par_in(3, 1:pix.nb_f) ;
                par_it = run_Opti(C, const_min, const_max, ...
                    option_opti, par_it) ;
                par_in(3, 1:pix.nb_f) = par_it ;
            otherwise
                error(['Fitting option ''', flag_fit, ...
                    ''' does not exist']) ;
        end
    end

    y_gauss = par_in(1, 1) ;
    x_gauss = par_in(2, 1) ;
    amp_gauss = par_in(3, :) ;
end