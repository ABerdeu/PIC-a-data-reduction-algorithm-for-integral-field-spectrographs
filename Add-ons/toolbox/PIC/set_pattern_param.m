%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to set the pattern parameters for a given wavelength (inversion
% operration of get_pattern_param)
% 
% Created: 08/31/2018 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #par_pat# the variable to update in array
% 
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
% #par_pat_out# updated cell array containing all the paremeters of the
%   pattern
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function par_pat_out = set_pattern_param(par_pat, par_pat_in, list_l, pat)

    nb_par = size(par_pat_in, 3) ;
    if nargin<4
        pat = 1 ;
    end
    par_pat_out = par_pat_in ;
    for par = 1:nb_par
        par_pat_out(pat, list_l, par) = par_pat(par,:) ;
    end
end