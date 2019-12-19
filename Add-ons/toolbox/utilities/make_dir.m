%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to create a directory if it is possible. If not it tries again
% or change its name
% 
% Created: 04/03/2017 (mm/dd/yyyy)
% Author: Anthony Berdeu (PhD student at CEA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #path_dir_in# the directory to create
%
%%%%%%%%
% Ouput
%%%%%%%%
% #path_dir_out# the actually created directory
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function path_dir_out = make_dir(path_dir_in)

    path_dir = path_dir_in ;
    try_count = 0 ; % To count the try before changing the name
    try_num = 1 ; % To count the number of tried names
    while try_count~=-1
        try
            % Trying to create the directory
            mkdir(path_dir) ;
            
            % Creation complete, exiting the loop
            path_dir_out = path_dir ;
            try_count=-1 ;
        catch exception
            try
                diary('./results/log_file.txt')
                disp(datestr(now,'yyyy_mm_dd_HH_MM_SS')) ;
                disp(exception) ;
                disp('Unable to create directory...') ;
                if try_count < 10
                    try_count = try_count + 1 ;
                    disp('Trying again...') ;
                else
                    try_count = 0 ;
                    path_dir = [path_dir_in, '_', num2str(try_num)] ;
                    try_num = try_num+1 ;
                    disp(['Changing name to :', path_dir]) ;
                end
                diary off
            catch
                disp(datestr(now,'yyyy_mm_dd_HH_MM_SS')) ;
                disp(exception) ;
                if try_count < 10
                    try_count = try_count + 1 ;
                    disp('Trying again...') ;
                else
                    try_count = 0 ;
                    path_dir = [path_dir_in, '_', num2str(try_num)] ;
                    try_num = try_num+1 ;
                    disp(['Changing name to :', path_dir]) ;
                end
            end
            pause(1) ;
        end
    end
end