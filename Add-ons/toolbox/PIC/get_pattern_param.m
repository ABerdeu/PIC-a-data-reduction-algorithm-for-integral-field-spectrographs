%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get the pattern parameters for a given wavelength
% 
% Created: 08/31/2018 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #par_pat_in# array containing all the paremeters of the pattern of size
%   [nb_pat, nb_l, nb_par]
%
% #list_l# list of lambda to extract
%
% #pat# pattern to extract (default 1: the first)
%
%%%%%%%%
% Ouput
%%%%%%%%
% #par_pat_out# the selected variable
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function par_pat_out = get_pattern_param(par_pat_in, list_l, pat)

    nb_par = size(par_pat_in, 3) ;
    nb_l = length(list_l) ;
    if nargin<3
        pat = 1 ;
    end
    par_pat_out = zeros(nb_par, nb_l) ;
    
    for par = 1:nb_par
        par_pat_out(par,:) = par_pat_in(pat, list_l, par) ;
    end
end