%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get the position of the local maximum of a region of interest
% in a picture
% 
% Created: 04/27/2018 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #pic# the picture
%
% #pos = [y, x]# the position around which the local maximum must be found
%
% #rad_ROI# radius of the region of interest to extract
%
% #flag_method# method to determine the local maximum
%   'max' -> pixel with the miximal value (default)
%   'centroid' -> pixel weitghed by their value
%
% #func# function to apply to the picture (default: identity)
%
% #flag_conv# several estimations until the position reaches convergence
% (default: false)
%
%%%%%%%%
% Ouput
%%%%%%%%
% #[m_y, m_x]# the pixel position of the local maximum
%
% #m# the maximum value
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [m_y, m_x, m] = find_local_maximum(pic, pos, rad_ROI, ...
    flag_method, func, flag_conv)

    % Initialization
    if ~isnumeric(pic)
        error('''find_local_maximum'' is only programmed for 2D arrays') ;
    end
    if size(pic,3) > 1
        error('''find_local_maximum'' is only programmed for 2D arrays') ;
    end
    
    [nb_y, nb_x] = size(pic) ;
    ind_y = (-rad_ROI:rad_ROI)+round(pos(1)) + floor(nb_y/2+1) ;
    ind_x = (-rad_ROI:rad_ROI)+round(pos(2)) + floor(nb_x/2+1) ;
    if nargin < 4
        flag_method = 'max' ;
    elseif isempty(flag_method)
        flag_method = 'max' ;
    end
        
    if nargin < 6
        flag_conv = false ;
    elseif isempty(flag_conv)
        flag_conv = false ;
    end
    
    % Insuring the extracted frame is in the picture
    ind_y = unique(max(1, min(nb_y, ind_y))) ;
    ind_x = unique(max(1, min(nb_x, ind_x))) ;
    
    % Maximum identification
    if ~isempty(ind_y) && ~isempty(ind_x)
        do = true ;
        m_y_old = Inf ;
        m_x_old = Inf ;
        
        while do
            % Extraction of the region of interest
            ROI_P = pic(ind_y, ind_x) ;

            % Application of func
            if nargin > 4
                if ~isempty(func)
                    ROI_P = func(ROI_P) ;
                end
            end
        
            % Determination of the local maximum
            switch flag_method
                case 'max' % Maximum pixel
                    m = max(ROI_P(:)) ;
                    [m_y, m_x] = ind2sub( ...
                        [length(ind_y), length(ind_x)], find(ROI_P==m)) ;
                    m_y = ind_y(m_y(1)) ;
                    m_x = ind_x(m_x(1));

                case 'centroid' % Centroid weighting
                    % Weighting
                    ROI_P = ROI_P/sum(ROI_P(:)) ;
                    [ind_x, ind_y] = meshgrid(ind_x, ind_y) ;

                    % Centroid 
                    m_y = ROI_P.*ind_y ;
                    m_y = sum(m_y(:)) ;
                    m_x = ROI_P.*ind_x ;
                    m_x = sum(m_x(:)) ;

                otherwise
                    error('Unknown method...') ;
            end
            
            % Next loop?
            do = flag_conv && ((m_y_old~=m_y)||(m_x_old~=m_x)) ;
            m_y_old = m_y ;
            m_x_old = m_x ;
            
            if do % Update the extracted positions
                ind_y = (-rad_ROI:rad_ROI)+round(m_y) ;
                ind_x = (-rad_ROI:rad_ROI)+round(m_x) ;

                % Insuring the extracted frame is in the picture
                ind_y = unique(max(1, min(nb_y, ind_y))) ;
                ind_x = unique(max(1, min(nb_x, ind_x))) ;
            end
        end
        
        % Going back to the frame coordinates
        m = pic(round(m_y), round(m_x)) ;
        m_y = m_y - floor(nb_y/2+1) ;
        m_x = m_x - floor(nb_x/2+1) ;
    else
        m = [] ;
        m_y = [] ;
        m_x = [] ;
    end
end