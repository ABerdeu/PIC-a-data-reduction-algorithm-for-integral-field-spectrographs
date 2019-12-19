%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to extract rows and columns in a matrix possibly stored in a bin
% file
% 
% Created: 03/18/2019 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #pic# picture in which the rows and columns must be extracted (can be a
% bin file)
%
% #list_ind = {list_i, list_j}# list of the index to extract
%
% #list_nb = [nb_y, nb_x]# Number of pixels in the global picture
%
%%%%%%%%
% Ouput
%%%%%%%%
% #pic_ij# extracted picture
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pic_ij = extract_ij(pic, list_ind, list_nb)

   if isnumeric(pic)
        % pic is a matrix
        pic_ij = pic(list_ind{1}, list_ind{2}) ;
    else
        % pic is a memmapfile
        pic_ij = extract_bin(pic, list_ind, list_nb) ;
    end
end

