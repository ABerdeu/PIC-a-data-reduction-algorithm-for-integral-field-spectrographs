%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get the bad pixels from a list of dark fitting
% 
% Created: 03/03/2019 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #[a_dark, b_dark]# maps of the coefficients fitted for the linear law
%
% #list_exp_time# list of the exposure times
%
% #list_avg# list of the averaged darks
%
% #list_var# list of the variance of the darks
%
% #tol_MAD# tolerance for bad pixels indentification
%
%%%%%%%%
% Ouput
%%%%%%%%
% #BP_mask# mask on the bad pixels
%
% #BP_mask_a# mask on the bad pixels identified with the slope statistic
%
% #BP_mask_res# mask on the bad pixels identified in the residues
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BP_mask, BP_mask_a, BP_mask_res] = ...
    get_BP_dark(a_dark, b_dark, list_exp_time, list_avg, list_var, tol_MAD)

    %% Initialization
    nb_dark = size(list_avg,3) ;

    %% Bad pixels on a_dark
    med_a = median(a_dark(:)) ;
    MAD_a = median(abs(a_dark(:)-med_a)) ;
    BP_mask_a = 1*(abs(a_dark-med_a)<tol_MAD*MAD_a) ;
    BP_mask_a = BP_mask_a.*(a_dark>0) ;

    % Percentage of bad pixels
    BP_percentage = sum(BP_mask_a(:))/numel(BP_mask_a(:)) ;

    %% Bad pixels on the cost
    BP_mask_res = ones(size(BP_mask_a)) ;
    
    % Weighted residues
    res = list_avg - ...
        a_dark.*reshape(list_exp_time, [1,1,nb_dark])-b_dark ; 
    
    % Cost at each pixel
    cost = sum(res(:,:,:).^2./list_var(:,:,:),3) ;
    
    % Removing already known bad pixels
    cost(~BP_mask_a) = 0 ;
    
    % Selection of the worst remaining pixels
    th_cost = quantile(cost(:), BP_percentage) ;
    BP_mask_res(cost>th_cost)=0 ;
    
    %% Global selection
    BP_mask = BP_mask_a.*BP_mask_res ;
end