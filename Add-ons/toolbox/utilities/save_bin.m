%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to save a binary file
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
% #var# the variable to save
%
% #file_format# the format to save
%
%%%%%%%%
% Ouput
%%%%%%%%
% #list_value# the list of the value of the matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save_bin(name_bin, var, file_format)
    count = 0 ;
    while count~=numel(var)
        fileID = fopen(name_bin,'w');
        count = fwrite(fileID, var, file_format);
        fclose(fileID);
        if count~=numel(var)
            try
                diary('./results/log_file.txt')
                disp(datestr(now,'yyyy_mm_dd_HH_MM_SS')) ;
                disp('Failure to save bin file... Trying again!') ;
                diary off
            catch
                disp(datestr(now,'yyyy_mm_dd_HH_MM_SS')) ;
                disp('Failure to save bin file... Trying again!') ;
            end
            pause(1) ;
        end
    end
end