%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get the bad pixels from a linear map
% 
% Created: 03/04/2019 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #[a_map, b_map]# maps of the coefficients fitted for the linear law
%
% #list_exp_time# list of the exposure times
%
% #list_avg# list of the averaged darks
%
% #list_var# list of the variance of the darks
%
% #tol_MAD# tolerance for bad pixels indentification
%
% #size_med_filter# size of the median filter to apply
%
%%%%%%%%
% Ouput
%%%%%%%%
% #BP_mask# mask on the bad pixels
%
% #BP_mask_a# mask on the bad pixels identified with the slope median
%
% #BP_mask_res# mask on the bad pixels identified in the residues
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BP_mask, BP_mask_a, BP_mask_res] = ...
    get_BP_mask(a_map, b_map, list_exp_time, list_avg, list_var, tol_MAD, ...
    size_med_filter)

    %% Initialization
    nb_acqui = size(list_avg,3) ;

    %% Bad pixels on a_map
    median_a = a_map - medfilt2(a_map, [size_med_filter,size_med_filter]) ;
    MAD_a = median(abs(median_a(:))) ;
    BP_mask_a = abs(median_a)<tol_MAD*MAD_a ;
    BP_mask_a = 1*BP_mask_a.*(a_map>0) ;

    % Percentage of bad pixels
    BP_percentage = sum(BP_mask_a(:))/numel(BP_mask_a(:)) ;

    %% Bad pixels on the cost
    BP_mask_res = ones(size(BP_mask_a)) ;
    
    % Weighted residues
    res = list_avg - ...
        a_map.*reshape(list_exp_time, [1,1,nb_acqui])-b_map ; 
    
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