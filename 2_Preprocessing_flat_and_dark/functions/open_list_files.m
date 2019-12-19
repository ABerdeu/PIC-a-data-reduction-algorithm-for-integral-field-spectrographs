%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to open a list of files with their exposure times
% 
% Created: 03/03/2019 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #file_path# path of the files
%
% #file_name# name of the files
%
%%%%%%%%
% Ouput
%%%%%%%%
% #nb_file# number of files
%
% #list_exp_time# list of the exposure times
%
% #list_avg# list of the averaged files
%
% #list_var# list of the variance of the files
%
% #list_file# list of the opened files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nb_file, list_exp_time, list_avg, list_var, list_file] = ...
    open_list_files(file_path, file_name)

    %% Initialization
    list_name = dir([file_path, '*', file_name, '*']) ;
    nb_file = length(list_name) ;
    return_list_file = nargout > 4 ;
    if return_list_file
        list_file = cell(nb_file, 1) ;
    end
    list_exp_time = zeros(nb_file,1) ;
    for file = 1:nb_file
        disp(['    Opening: ', num2str(file), '/', num2str(nb_file)]) ;
        
        % Name of the file
        name_f = [file_path, list_name(file).name] ;
        % Opening data
        file_f = fitsread(name_f) ;
        if return_list_file
            list_file{file} = file_f ;
        end
        
        % Exposure time
        list_exp_time(file) = get_fits_exposure(name_f) ;
        
        
        % Initializing list_avg and list_var
        if file==1
            list_avg = zeros(size(file_f,1), ...
                size(file_f,2), ...
                nb_file) ;
            list_var = list_avg ;
        end
        
        % Average
        list_avg(:,:,file) = mean(file_f, 3) ;
        
        % Variance
        list_var(:,:,file) = var(file_f, [], 3) ;
    end
end