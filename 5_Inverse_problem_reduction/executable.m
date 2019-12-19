%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This macro is a procedure to run an inverse problem approach on IFS data
% knowing the projection model
%
% Created: 06/01/2018 (mm/dd/yyyy)
% Modified: 08/08/2018 (mm/dd/yyyy) Taking into account the transmission
% given by the flat. Reweighted L2 norm
% Modified: 09/07/2018 (mm/dd/yyyy) Introducing spectral response and
% transmission calibrated with the flat
% Modified: 10/16/2018 (mm/dd/yyyy) Taking into account the full field flat
% and the Poisson noise in the robust penalization
% Modified: 02/05/2019 (mm/dd/yyyy) Taking into account the background
% calibration
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
% Modified: 10/14/2019 (mm/dd/yyyy) On-line deposit
% Author: Anthony Berdeu (Centre de Recherche en Astrophysique de Lyon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all ;
close all ;
clc ;

% Absolute path
abs_path = pwd ;
abs_path = [abs_path, '\'] ;

% Load functions and packages
restoredefaultpath
addpath(genpath('./functions')) ;
addpath('./procedures') ;

% Loading toolboxes
addpath(genpath('../Add-ons/')) ;


%% Set parameters
run parameters
save_path = [save_path, ...
    datestr(now,'yyyy_mm_dd_HH_MM_SS'),'/'] ;
save_path = make_dir(save_path) ;

make_dir([save_path, '/IFS_sim']) ;
make_dir([save_path, '/Obj']) ;
make_dir([save_path, '/Res']) ;
make_dir([save_path, '/W_coef']) ;
make_dir([save_path, '/IFS_BP']) ;

% Saving the running codes
save_path_aux = make_dir([save_path, 'code/']) ;
copyfile('./executable.m', save_path_aux) ;
save_path_aux = make_dir([save_path, 'code/procedures/']) ;
copyfile('./procedures',save_path_aux) ;


%% Initialization
global IFS_BP_W
run initialization

%% Iterative reconstruction
% Initialization with the transposition of the direct model
% hyp_cub = 0*Forward_Mod'*(IFS_data-IFS_BG) ;
% Initialization with 0
hyp_cub = 0*Forward_Mod'*(IFS_data-IFS_BG) ;
IFS_sim = Forward_Mod*hyp_cub+IFS_BG ;
Res = IFS_data-IFS_sim ;
save_fits(hyp_cub, 'Obj/Obj_it_0', save_path) ;
save_fits(IFS_sim, 'IFS_sim/IFS_sim_it_0', save_path) ;
save_fits(Res, 'Res/Res_it_0', save_path) ;

for it = 1:nb_res_estimation
    disp(['Global iteration: ', num2str(it), '/', ...
        num2str(nb_res_estimation)]) ;
    
    %% Building cost function
    % Data fidelity
    if it == 1
        option_robust = [] ;
        option_robust.method = option_opti.RL2_method ;
        
        % Rough suppression of the bad pixel by a median filter
        var_s = medfilt2(IFS_data, med_filter_IFS) ;
        save_fits(var_s, 'Corrected_IFS_data_medfilter', save_path) ;
        var_s = max(0, var_s) ;
        var_s = eta*var_s+var_0 ;
        
        % Replacing negative values
        if sum(var_s(:)<=0)
            warning(['There are ', num2str(sum(var_s(:)<=0)), ...
                ' elements of the variance lower than 0...']) ;
            warning(['There are set to ', num2str(eps)]) ;
            var_s(var_s<=0) = eps ;
        end
        option_robust.flag_s = IFS_BP_W.* ...
            list_scaling_factor(it).*(var_s).^0.5 ;
        option_robust.noise_model = 'none' ;
        Data_fid = CostRobustPenalization(Forward_Mod, ...
            IFS_data-IFS_BG, ...
            option_robust) ;
        
        % Saving initial weights
        W_coef = Data_fid.computeW_(hyp_cub) ;
        save_fits(W_coef, 'W_coef/W_coef_it_0', save_path) ;
    else
        option_robust = [] ;
        option_robust.method = option_opti.RL2_method ;
        if list_Poisson_flag(it) % The Poisson noise is dynamic
            option_robust.flag_s = IFS_BP_W.*list_scaling_factor(it) ;
            option_robust.var_0 = var_0+eta*IFS_BG ;
            option_robust.eta = eta ;
            option_robust.noise_model = 'Poisson' ;
        else % The Poisson noise is in the scaling of the residues
            option_robust.flag_s = IFS_BP_W.*list_scaling_factor(it).* ...
                (eta.*IFS_sim+var_0).^0.5 ;
            option_robust.noise_model = 'none' ;
        end
        Data_fid = CostRobustPenalization(Forward_Mod, ...
            IFS_data-IFS_BG, ...
            option_robust) ;
    end
    
    % Global cost function
    Cost = 1/prod([pix_IFS.nb_y, pix_IFS.nb_x])*Data_fid ;
    if list_mu_TV(it)>0
        Cost = Cost + ...
            list_mu_TV(it)/list_scaling_factor(it)^2/prod( ...
            [pix_hypcube.nb_y, pix_hypcube.nb_x, pix_hypcube.nb_l] ...
            )*Cost_TV ;
    end
    if list_L2_grad_lambda(it)>0
        Cost = Cost + ...
            list_L2_grad_lambda(it)/list_scaling_factor(it)^2/prod( ...
            [pix_hypcube.nb_y, pix_hypcube.nb_x, pix_hypcube.nb_l] ...
            )*Cost_L2_grad_lambda ;
    end
    
    %% Optimization
    hyp_cub = run_Opti(Cost, hyp_cub_min, hyp_cub_max, option_opti, ...
        hyp_cub) ;
    
    %% Saving iteration
    % Simulation and residues
    IFS_sim = Forward_Mod*hyp_cub ;
    Res = IFS_data-IFS_sim-IFS_BG ;
    W_coef = Data_fid.computeW_(hyp_cub) ;
    
    % Saving
    save_fits(IFS_sim, ['IFS_sim/IFS_sim_it_', num2str(it)], save_path) ;
    save_fits(W_coef, ['W_coef/W_coef_it_', num2str(it)], ...
        save_path) ;
    save_fits(Res, ['Res/Res_it_', num2str(it)], save_path) ;

    % Reconstructed object
    save_fits(hyp_cub, ['Obj/Obj_it_', num2str(it)], save_path) ;
    
    % Adding background in the model
    IFS_sim = IFS_sim+IFS_BG ;
    
    %% Rough suppression of the bad pixel by a median filter
    if list_median_flag(it)
        hyp_cub = medfilt3(hyp_cub, med_filter_cube) ;
        
        % Simulation and residues
        IFS_sim = Forward_Mod*hyp_cub ;
        Res = IFS_data-IFS_sim-IFS_BG ;
        
        % Saving
        save_fits(hyp_cub, ['Obj/Obj_it_', num2str(it), ...
            '_medfilter'], save_path) ;
        W_coef = Data_fid.computeW_(hyp_cub) ;
        save_fits(W_coef, ['W_coef/W_coef_it_', num2str(it), ...
            '_medfilter'], save_path) ;
        save_fits(IFS_sim, ['IFS_sim/IFS_sim_it_', num2str(it), ...
            '_medfilter'], save_path) ;
        save_fits(Res, ['Res/Res_it_', num2str(it), ...
            '_medfilter'], save_path) ;
        flag_BP = 1./IFS_BP_W ;
        flag_BP(IFS_BP==0) = 1 ;
        save_fits(uint8(flag_BP), ...
            ['IFS_BP/IFS_BP_it_', num2str(it), ...
            '_medfilter'], save_path) ;
        
        % Adding background in the model
        IFS_sim = IFS_sim+IFS_BG ;
    end
    
    %% Thresholding the weighting coefficients
    if list_threshold_W_coef(it)
        % Computing weighting coeffcients without the prior on the bad
        % pixels
        option_robust = [] ;
        option_robust.method = option_opti.RL2_method ;
        option_robust.flag_s = list_scaling_factor(it).* ...
            (eta.*IFS_sim+var_0).^0.5 ;
        option_robust.noise_model = 'none' ;
        Data_fid = CostRobustPenalization(Forward_Mod, ...
            IFS_data-IFS_BG, ...
            option_robust) ;
        W_coef = Data_fid.computeW_(hyp_cub) ;
        
        % New prior on the bad pixels
        IFS_BP_W = ones(pix_IFS.nb_y, pix_IFS.nb_x) ;
        IFS_BP_W(W_coef<th_W_coef) = Inf ;
        IFS_BP_W(IFS_BP==0) = Inf ;
    end
    flag_BP = 1./IFS_BP_W ;
    flag_BP(IFS_BP==0) = 1 ;
    save_fits(uint8(flag_BP), ...
        ['IFS_BP/IFS_BP_it_', num2str(it)], save_path) ;
end