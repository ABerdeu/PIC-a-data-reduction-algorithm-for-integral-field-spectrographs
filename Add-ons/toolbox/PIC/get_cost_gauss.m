%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get the cost function to determine the parameters of a
% gaussian pattern with some known
% parameters
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
% #list_y# list of the positions on the y-axis on which the gaussian
% pattern is extracted
% 
% #list_x# list of the position on the x-axis on which the gaussian
% pattern is extracted
% 
% #flag_s# method to scale the residue before applying the
%   objective function (default: none)
%   ->  'none'	-> no scaling applied
%   ->  'MAD'   -> median absolute deviation
%   ->  numeric -> scaling
%
% #RL2_method# method for the reweighted least square
%
% #y_c# coordinate of the center of the gaussian pattern on the
% y-axis (if empty, it becomes a variable to estmate)
% 
% #x_c# coordinate of the center of the gaussian pattern on the
% x-axis (if empty, it becomes a variable to estmate)
% 
% #theta# orientation (in degree) of the gaussian pattern
% (if empty, it becomes a variable to estmate)
% 
% #sig_par# elongation of the gaussian pattern along the
% direction given by theta (if empty, it becomes a variable to give as
% an input)
% 
% #sig_perp# elongation of the gaussian pattern along the
% perpendicular direction given by theta (if empty, it becomes a
% variable to estmate)
% 
% #amp# amplitude of the gaussian pattern (if empty, it becomes a
% variable to estmate)
% 
% #offset# offset of the gaussian pattern (if empty, it becomes a
% variable to estmate)
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
function [C] = get_cost_gauss(pic, list_y, list_x, flag_s, ...
    RL2_method, y_c, x_c, theta, sig_par, sig_perp, amp, offset, list_bool)


    %% Initialization
    [nb_y, nb_x, nb_frame] = size(pic) ;
    if nargin < 13
        list_bool = true ;
    end
    if sum(list_bool) == length(list_bool) % All the variable are the same
            % in all frame, no selection needed
        list_bool = [] ;
    end
    nb_var = length(list_bool) ;
    
    % Checking if inputs are depending of the frame
    y_fd = length(y_c)>1 ;
    x_fd = length(x_c)>1 ;
    theta_fd = length(theta)>1 ;
    sig_par_fd = length(sig_par)>1 ;
    sig_perp_fd = length(sig_perp)>1 ;
    amp_fd = length(amp)>1 ;
    offset_fd = length(offset)>1 ;


    %% Building cost function
    for f = 1:nb_frame
        % Variable values in the frame
        y_c_f = selection_variable(y_c, f, y_fd) ;
        x_c_f = selection_variable(x_c, f, x_fd) ;
        theta_f = selection_variable(theta, f, theta_fd) ;
        sig_par_f = selection_variable(sig_par, f, sig_par_fd) ;
        sig_perp_f = selection_variable(sig_perp, f, sig_perp_fd) ;
        amp_f = selection_variable(amp, f, amp_fd) ;
        offset_f = selection_variable(offset, f, offset_fd) ;
        
        % Extracting picture
        pic_f = pic(1:nb_y, 1:nb_x, f) ;

        % Computing the scaling factor
        s = get_scaling_factor(pic_f, flag_s) ;

        % Gaussian operator
        Gauss = OpGauss(list_y, list_x, ...
            y_c_f, x_c_f, theta_f, sig_par_f, sig_perp_f, amp_f, ...
            offset_f) ;

        % Cost for the frame f
        if isempty(list_bool)
            C_aux = CostReweightedL2(Gauss, pic_f, RL2_method, ...
                [], s) ;
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
            C_aux = CostReweightedL2(Gauss*select_var, ...
                pic_f, RL2_method, [], s) ;
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