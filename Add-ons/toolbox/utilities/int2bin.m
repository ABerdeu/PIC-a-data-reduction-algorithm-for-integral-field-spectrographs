%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to convert a positive integer to a binary number
%
% Created: 04/06/2018 (mm/dd/yyyy)
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #m# a positive integer
%
% #n_in# the wanted length of the binary vector (optional)
%
%%%%%%%%
% Ouput
%%%%%%%%
% #b# m written in a binary mode
%
% #n_out# the output length of the binary vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [b, n_out] = int2bin(m, n_in)

    %% Initialization
    n_out = floor(log(m)/log(2))+1 ;
    if nargin>1
        n_out = max(n_out, n_in) ;
    end

    b = zeros(n_out,1) ;
    r = m ;
    
    %% Loop on the powers of 2
    for p = n_out-1:-1:0
        q = floor(r/2^p) ;
        r = r-q*2^p ;
        b(p+1) = q ;
    end
end