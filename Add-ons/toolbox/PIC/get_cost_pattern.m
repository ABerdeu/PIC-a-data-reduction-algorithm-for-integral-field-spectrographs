%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get the cost function to determine the parameters of a
% specific pattern with some known parameters
% 
% Created: 05/02/2018 (mm/dd/yyyy)
% Modified: 08/31/2018 (mm/dd/yyyy) Inspired from get_cost_gaussian to
% shift to different patterns
% Modified: 11/06/2018 (mm/dd/yyyy) New definition of the robust
% penalization
% Modfied: 02/26/2019 (mm/dd/yyyy) Hexagon pattern
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #pic# picture in which the pattern is extracted
%
% #list_y# list of the positions on the y-axis on which the gaussian
% pattern is extracted
% 
% #list_x# list of the position on the x-axis on which the gaussian
% pattern is extracted
% 
% #pattern_model# flags to specify the pattern to fit
%   pattern_model.flag_norm -> Normalized?
%   pattern_model.flag_profile -> Gaussian / Hexagon / Moffat
%   pattern_model.oversampling -> oversampling of the pattern (default: 1)
%
% #pattern_car# cell containing the list of the pattern characteristic 
% (positions, amplitude, offset...) (if a field is empty, it becomes a
% variable to estimate)
% 
% #option_robust# options for the robust penalization
% 
% #list_bool# list of booleans; for eaxh variable to estimate it precises
% if it the same for all frame (true) or different (false). Default: true.
%
%%%%%%%%
% Ouput
%%%%%%%%
% #C# Global cost
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C] = get_cost_pattern(pic, list_y, list_x, pattern_model, ...
    pattern_car, option_robust, list_bool)

    %% Initialization
    [nb_y, nb_x, nb_frame] = size(pic) ;
    if nargin < 7
        list_bool = true ;
    end
    if sum(list_bool) == length(list_bool) % All the variable are the same
            % in all frame, no selection needed
        list_bool = [] ;
    end
    % Number of variables which are the same in all the frames
    nb_var = length(list_bool) ;
    
    % Number of inputs to describe the pattern
    nb_input = length(pattern_car) ;
    
    % Checking if inputs are depending of the frame
    input_fd = false(nb_input, 1) ;
    for in = 1:nb_input
        input_fd(in) = length(pattern_car{in})>1 ;
    end

    %% Building cost function
    for f = 1:nb_frame
        % Variable values in the frame
        input_f = cell(nb_input, 1) ;
        for in = 1:nb_input
            input_f{in} = selection_variable(pattern_car{in}, f, ...
                input_fd(in)) ;
        end
        
        % Extracting picture
        pic_f = pic(1:nb_y, 1:nb_x, f) ;

        % Computing the scaling factor
        option_robust.flag_s = ...
            get_scaling_residues_factor(pic_f, option_robust.flag_s) ;
        
        if isfield(pattern_model, 'oversampling') && ...
                pattern_model.oversampling>1
            [OpDS, list_y_OS, list_x_OS] = get_downsampling_operator( ...
                pattern_model.oversampling, list_y, list_x) ;
        else
            OpDS = 1 ;
            list_y_OS = list_y ;
            list_x_OS = list_x ;
        end

        % Getting the pattern operator
        switch pattern_model.flag_profile
            case 'Gaussian'
                % Input: y_c, x_c, sig, amp, offset
                Pattern = OpDS*OpGauss(list_y_OS, list_x_OS, ...
                    input_f{1}, input_f{2}, input_f{3}, input_f{4}, ...
                    input_f{5}, pattern_model.flag_norm) ;
            case 'Hexagon'
                % Input: y_c, x_c, sig, amp, offset
                Pattern = OpDS*OpPSF_Hex(list_y_OS, list_x_OS, ...
                    input_f{1}, input_f{2}, pattern_model.theta, ...
                    input_f{3}, input_f{4}, ...
                    input_f{5}, pattern_model.flag_norm) ;
            case 'Moffat'
                % Input: y_c, x_c, alpha, beta, amp, offset
                Pattern = OpDS*OpMoffat(list_y_OS, list_x_OS, ...
                    input_f{1}, input_f{2}, input_f{3}, input_f{4}, ...
                    input_f{5}, input_f{6}, pattern_model.flag_norm) ;
            otherwise
                error([pattern_model.flag_profile, ...
                    ': Unkown model for the pattern...']) ;
        end

        % Cost for the frame f
        if isempty(list_bool)
            C_aux = CostRobustPenalization(Pattern, pic_f, option_robust) ;
        else
            sel = false(nb_var, nb_frame) ;
            for var = 1:nb_var
                if list_bool(var) % The variable is the same in all frame
                    sel(var,1) = true ;
                else
                    sel(var,f) = true ;
                end
            end
            select_var = LinOpSelector(sel) ;
            C_aux = CostRobustPenalization(Pattern*select_var, ...
                pic_f, option_robust) ;
        end
        
        % Global cost
        if f==1
            C = C_aux ;
        else
            C = C + C_aux ;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to select the correct variable for a given frame
% 
% Created: 07/04/2018 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #var_in# the list of input variable
%
% #ind# the selected index (if possible)
%
% #frame_dep# is the variable depending of the frame?
%
%%%%%%%%
% Ouput
%%%%%%%%
% #sel_var# the selected variable
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sel_var = selection_variable(var_in, ind, frame_dep)
    if frame_dep
        sel_var = var_in(ind) ;
    else
        sel_var = var_in ;
    end
end