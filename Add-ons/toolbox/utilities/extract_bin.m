%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to extract a part of a binary file
% 
% Created: 07/04/2018 (mm/dd/yyyy)
% Author: Anthony Berdeu (PhD student at CEA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #memmapfile# pointer to the matrix
%
% #list_ind# the list of index to extract
%
% #list_dim# the list of the dimensions of the matrix stored in memmapfile
%
%%%%%%%%
% Ouput
%%%%%%%%
% #sub_part# the extracted subpart
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sub_part] = extract_bin(memmapfile, list_ind, list_dim)
    %% Initialization
    nb_dim = length(list_dim) ;
    dim_out = zeros(1, nb_dim) ;
    for dim = 1:nb_dim
        dim_out(dim) = length(list_ind{dim}) ;
    end
    
    %% List of indexes to extract
    ind_ext = list_ind{nb_dim} ;
    for dim = (nb_dim-1):-1:1
        [ind_1, ind_2] = meshgrid(ind_ext, list_ind{dim}) ;
        ind_ext = list_dim(dim)*(ind_1(:)-1) +  ind_2(:) ;
    end

    %% Extraction
    sub_part = memmapfile.Data(ind_ext) ;
    sub_part = reshape(sub_part, dim_out) ;
end