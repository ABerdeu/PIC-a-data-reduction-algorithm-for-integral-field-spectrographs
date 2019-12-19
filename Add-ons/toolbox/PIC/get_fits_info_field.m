%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get the information present in a given field of a fits info
% files
% 
% Created: 01/24/2019 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #fits_info# header of the fits file
%
% #target_field# String of character of the targeted field
% 
% #delimiter# delimited between which the information is contained
% 
%%%%%%%%
% Ouput
%%%%%%%%
% #field# The wanted field
%
% #field_lign# The wanted field lign
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [field, field_lign] = get_fits_info_field(fits_info, ...
    target_field, delimiter)
    
    % List of the fields
    list_field = fits_info.PrimaryData.Keywords ;
    
    % The fields are in the third column
    list_field = list_field(:,3) ;
    nb_field = length(list_field) ;
    
    % Finding field
    lign = 1 ;
    found = false ;
    while lign < (nb_field+1) && ~found
        pos_field = strfind(list_field{lign}, target_field) ;
        
        if isempty(pos_field)
            lign = lign+1 ;
        else
            found = true ;
        end
    end
    
    % Extracting field
    field_lign = list_field{lign}((pos_field+length(target_field)):end) ;
    field = field_lign ;
    field = strsplit(field,delimiter) ;
    field = field{2} ;
end