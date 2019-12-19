%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to read a full binary file
% 
% Created: 03/08/2017 (mm/dd/yyyy)
% Modified: 04/11/2017 (mm/dd/yyyy) Taking in account the format
% Author: Anthony Berdeu (PhD student at CEA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #name_bin# the name of the binary file (full path needed)
%
% #var# the variable to read
%
% #file_format# the format to read
%
%%%%%%%%
% Ouput
%%%%%%%%
% #val_bin# the list of the value of the matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [val_bin] = read_bin(name_bin, file_format)
    fileID = fopen(name_bin,'r');
    val_bin = fread(fileID, file_format);
    fclose(fileID);
end