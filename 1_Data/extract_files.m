%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This macro is a procedure to extract the files and renaming them
%
% Created: 01/24/2019 (mm/dd/yyyy)
% Author: Anthony Berdeu (Laboratoire Hubert Curien)
% Modified: 10/14/2019 (mm/dd/yyyy) On-line deposit
% Author: Anthony Berdeu (Centre de Recherche en Astrophysique de Lyon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all ;
close all ;
clc ;

% Absolute path
abs_path = pwd ;
abs_path = [abs_path, '\'] ;

% Load functions and packages
restoredefaultpath

% Loading toolboxes
addpath(genpath('../Add-ons/')) ;

%% Parameters
% path to the raw data downloaded from the data center
folder_in = 'path_raw_data_in\' ;

% path to the output folder to store the renamed files
folder_out = 'path_renamed_data_out\' ;

%% Initialization
% Output folder
if ~exist(folder_out, 'dir')
    make_dir(folder_out) ;
end

% List of fits files
list_fits = dir([folder_in, '*.fits']) ;
nb_fits = length(list_fits) ;

%% Renaming files
for file_i = 1:nb_fits
    disp(['Renaming file ', num2str(file_i), '/', num2str(nb_fits)]) ;
    
    % Suffix
    suffix = get_fits_suffix([folder_in, list_fits(file_i).name]) ;
    
    % New name
    [~, name_out, ~] = fileparts([folder_in, list_fits(file_i).name]) ;
    
    % Copying
    copyfile([folder_in, list_fits(file_i).name], ...
        [folder_out, [name_out, suffix, '.fits']]) ;
end