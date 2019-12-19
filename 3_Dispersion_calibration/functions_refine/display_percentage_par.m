%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to display a percentage in a parallel loop. It assumes that it
% has been initialized with a first display of the percentage 000.00 %
% 
% Created: 03/20/2019 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #path_count# path to save the count
%
% #delta_count# every delta_count, the percentage is displayed
%
% #nb_tot_count# number total of counts
%
%%%%%%%%
% Ouput
%%%%%%%%
% #count# the count
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function count = display_percentage_par(path_count, ...
    delta_count, nb_tot_count)

    %% Opening the counts
    failed = 1 ;
    while failed
        try
            count = fitsread([path_count, 'count.fits']) ;
            % Update of the counts
            save_fits(count+1, 'count', path_count) ;
            failed = 0 ;
        catch
            count = 1 ;
            failed = 1 ;
        end
    end
    
    %% Display the percetage if needed
    if mod(count, delta_count) == 0 % Display percentage?
        fprintf(repmat('\b', 1, 11));
        disp([': ', num2str(count/nb_tot_count*100, '%06.2f'), ' %']) ;
    end
end

