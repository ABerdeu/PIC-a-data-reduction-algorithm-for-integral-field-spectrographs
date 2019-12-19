%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This macro is a procedure to analyze the centroid of a coronographic
% acquisition
%
% Created: 10/30/2018 (mm/dd/yyyy)
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


%% Set parameters
% Path to the data and name of the file
data.path = './results/path_to_the_data/' ;
data.name = 'data_name.fits' ;

% First and last frames to find the star spots
l_first = 2 ;
l_end = 44 ;

% Lambda
list_lambda = ['../4_PIC_operators_declaration/', ...
    'results/path_to_the_model/list_lambda.txt'] ;

% Guess by the user of the gaussian parameter
min_size = 1 ;
guess_size = 2 ;
flag_norm = true ; % Is the pattern normalized?

% nb_cal_pos = 10 ; % < nb_cal_par
nb_cal = 10 ;
tol_center = 0.5 ; % Tolerance for the center positionning
nb_tol_plot = 5 ; % Number of tolerance to plot
flag_offset = false ; % Offset in the linear fit of the wavelength scaling

list_cal_pat = {'amp', 'pos'} ; % Calibration technique for the
    % pattern fitting

% Option of the optimizer for the position and amplitude
paralleling = true ;


% Options for the optimization of the positions and size
option_opti_amp_pos.RL2_method = 'Cauchy' ;
option_opti_amp_pos.noise_model = 'none' ;
option_opti_amp_pos.eta = 0 ;
option_opti_amp_pos.var_0 = 1 ;
option_opti_amp_pos.flag_s = 100000 ;

% Option of the optimizer for the position and amplitude
option_opti_amp_pos.method = 'fminsearch' ;
option_opti_amp_pos.verbose = false ;
option_opti_amp_pos.maxiter = 250 ; % Maximal number of iterations for the
    % optimizer
option_opti_amp_pos.memoizeOpts = true ; % Memorization of the costs
    % computations to save time in the iteration process

% Option of the optimizer for the pattern parameters in the initialization
option_opti_par = option_opti_amp_pos ;
option_opti_par.maxiter = 250 ;
option_opti_par.verbose = false ;
option_opti_par.method = 'fminsearch' ;

% Options for the square fitting
option_opti.method = 'VMLMB' ;
option_opti.verbose = false ;
option_opti.maxiter = 50 ; % Maximal number of iterations for the
    % optimizer
option_opti.memoizeOpts = true ; % Memorization of the costs computations
    % to save time in the iteration process

% Saving png
flag_png = true ;

% Display option
color_scale = [0, 25000] ;
nb_color = 1024 ;
color_map = jet(nb_color) ;
fig_pos = [250, 100] ;
fig_size = 900*[1.5, 1] ;
rad_ROI = 6 ;

% Marker and font size
MarkerSize = 10 ;
LineWidth = 2.5 ;
FontSize = 20 ;
FontSize_axis = 15 ;
scale_plot = 2.5 ;


%% Opening data
save_path = [data.path, 'centroid_analysis/'] ;
if ~exist(save_path, 'dir')
    save_path = make_dir(save_path) ;
end
if ~exist([save_path, 'square_fit/'], 'dir')
    make_dir([save_path, 'square_fit/']) ;
end
copyfile('./analyze_centroid.m', save_path) ;

% Opening data cube
hyp_cube = fitsread([data.path, data.name]) ;
hyp_cube_sim = hyp_cube*0 ;


%% Initialization
% Gaussian pattern
pattern_model.flag_norm = true ; % Normalized?
pattern_model.oversampling = 1 ; % Scale to oversample the pattern to fit 
    % (must be an integer)
pattern_model.flag_profile = 'Gaussian' ;


[pix_cube.nb_y, pix_cube.nb_x, pix_cube.nb_l] = size(hyp_cube) ;
pix_cube.dx = 1 ;
pix_cube.dy = 1 ;
pix_cube.nb_f = 1 ;
list_lambda = load(list_lambda, '-ascii') ;

% Corner
if exist([save_path, 'list_y.mat'], 'file')
    hyp_cube_sim = fitsread([save_path, 'hyp_cube_sim.fits']) ;
    load([save_path, 'list_y.mat']) ;
    load([save_path, 'list_x.mat']) ;
    load([save_path, 'list_sig.mat']) ;
    load([save_path, 'list_amp.mat']) ;
    load([save_path, 'list_off.mat']) ;
    load_corner = true ;
else
    list_y = zeros(pix_cube.nb_l, 4) ;
    list_x = zeros(pix_cube.nb_l, 4) ;
    list_amp = zeros(pix_cube.nb_l, 4) ;
    list_sig = guess_size*ones(pix_cube.nb_l, 4) ;
    list_off = zeros(pix_cube.nb_l, 1) ;
    load_corner = false ;
end

% Centroid
if exist([save_path, 'list_centroid.mat'], 'file')
    load([save_path, 'list_centroid.mat']) ;
    load([save_path, 'list_diag.mat']) ;
    load([save_path, 'list_theta.mat']) ;
    load_centroid = true ;
else
    list_centroid = zeros(pix_cube.nb_l, 2) ;
    list_diag = zeros(pix_cube.nb_l, 1) ;
    list_theta = zeros(pix_cube.nb_l, 1) ;
    load_centroid = false ;
end

%% Opening parralel workers
if paralleling && ~load_corner
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        try
            poolobj = parpool ;
            disp(['Connected to ', num2str(poolobj.NumWorkers), ...
            ' workers.']) ;
        catch
            disp('Unable to start workers, license non available...') ;
            disp('No paralleling jobs...') ;
        end
    else
        disp(['Already connected to ', num2str(poolobj.NumWorkers), ...
            ' workers.']) ;
    end
end


%% First frame analysis
if ~ load_corner
    first_frame = true ;
    fig_1 = figure(1) ;
    l = l_first ;
    frame_l = hyp_cube(:,:,l) ;
    med_frame = median(frame_l(:)) ;
    list_off(l) = med_frame ;
    while first_frame
        %% Plot figures
        imshow(frame_l, color_scale) ;
        colormap(gca, color_map) ;
        title('Select two opposite corners...', 'FontSize', FontSize) ;
        set(findall(gca,'-property','Interpreter'),'Interpreter','Latex') ;
        set(findall(gcf,'-property','Interpreter'),'Interpreter','Latex') ;
        set(fig_1,'rend','painters','pos', [fig_pos, fig_size]) ;

        %% Select first 2 points
        for corner = 1:2:3
            % User pointing
            [list_x(l, corner), list_y(l, corner)] = ...
                my_ginput(1, [0,0,0]) ;
            ind_y = round(list_y(l, corner)) ;
            ind_x = round(list_x(l, corner)) ;
            list_amp(l, corner) = frame_l(ind_y, ind_x) ;

            list_y(l, corner) = list_y(l, corner) - ...
                floor(pix_cube.nb_y/2+1) ;
            list_x(l, corner) = list_x(l, corner) - ...
                floor(pix_cube.nb_x/2+1) ;

            % Display result
            plot_rectangle(fig_1, ...
                [list_y(l, corner), list_x(l, corner)], ...
                scale_plot*list_sig(l, corner)*[1, 1], 0, pix_cube, ...
                [1,1,1]) ;

            % Set marker size
            list_plot = get(gca(fig_1),'Children') ;
            for h = 1:length(list_plot)
                if strcmp('line', get(list_plot(h), 'Type'))
                    set(list_plot(h), 'MarkerSize', MarkerSize, ...
                        'LineWidth', LineWidth) ;
                end
            end
            pause(0.1) ;
        end

        %% Guessing first square
        list_centroid(l, :) = [sum(list_y(l, [1,3]))/2 , ...
            sum(list_x(l, [1,3]))/2] ;
        list_diag(l) = 0.5 * ...
            ((list_y(l,3)-list_y(l,1))^2+(list_x(l,3)-list_x(l,1))^2)^0.5 ;
        list_theta(l) = mod( ...
            atand((list_y(l,3)-list_y(l,1))/(list_x(l,3)-list_x(l,1))), ...
            90) ;
        list_corner = list_centroid(l, :) + ...
            list_diag(l)*[sind(list_theta(l)+(0:90:270)'), ...
            cosd(list_theta(l)+(0:90:270)')] ;
        list_y(l, :) = list_corner(:,1) ;
        list_x(l, :) = list_corner(:,2) ;

        % Display
        clf(fig_1) ;
        imshow(frame_l, color_scale) ;
        colormap(gca, color_map) ;
        title('Refining the first guess...', 'FontSize', FontSize) ;
        set(findall(gca,'-property','Interpreter'),'Interpreter','Latex') ;
        set(findall(gcf,'-property','Interpreter'),'Interpreter','Latex') ;
        plot_rectangle(fig_1, ...
            [list_y(l, :)', list_x(l, :)'], ...
            scale_plot*list_sig(l, :)'.*[1, 1], 0, pix_cube, [1,1,1]) ;
        % Set marker size
        list_plot = get(gca(fig_1),'Children') ;
        for h = 1:length(list_plot)
            if strcmp('line', get(list_plot(h), 'Type'))
                set(list_plot(h), 'MarkerSize', MarkerSize, ...
                    'LineWidth', LineWidth) ;
            end
        end
        set(fig_1,'rend','painters','pos', [fig_pos, fig_size]) ;
        pause(0.1) ;

        %% Refining the first square
        % Loop on the corners
        if paralleling
            for corner = 1:4
                for cal = 1:nb_cal
                    % Refining position and amplitude
                    for cc = 1:length(list_cal_pat)
                        [list_y(l, corner), list_x(l, corner), ...
                            list_amp(l, corner)] = ...
                            find_pattern_pos_amp(frame_l, ...
                            [list_y(l, corner); list_x(l, corner); ...
                            list_amp(l, corner)], ...
                            rad_ROI, pix_cube, option_opti_amp_pos, ...
                            pattern_model, ...
                            list_sig(l, corner), med_frame, ...
                            list_cal_pat{cc}) ;
                    end

                    % Refining pattern parameter
                    list_sig(l, corner) = find_pattern_par(frame_l, ...
                        list_sig(l, corner), ...
                        rad_ROI, pix_cube, option_opti_par, ...
                        pattern_model, ...
                        list_y(l, corner), list_x(l, corner), ...
                        list_amp(l, corner), med_frame) ;
                end
            end
        else
            for corner = 1:4
                for cal = 1:nb_cal
                    % Refining position and amplitude
                    for cc = 1:length(list_cal_pat)
                        [list_y(l, corner), list_x(l, corner), ...
                            list_amp(l, corner)] = ...
                            find_pattern_pos_amp(frame_l, ...
                            [list_y(l, corner); list_x(l, corner); ...
                            list_amp(l, corner)], ...
                            rad_ROI, pix_cube, option_opti_amp_pos, ...
                            pattern_model, ...
                            list_sig(l, corner), med_frame, ...
                            list_cal_pat{cc}) ;
                    end

                    % Refining pattern parameter
                    list_sig(l, corner) = find_pattern_par(frame_l, ...
                        list_sig(l, corner), ...
                        rad_ROI, pix_cube, option_opti_par, ...
                        pattern_model, ...
                        list_y(l, corner), list_x(l, corner), ...
                        list_amp(l, corner), med_frame) ;
                end
            end
        end

        % Display
        clf(fig_1) ;
        imshow(frame_l, color_scale) ;
        colormap(gca, color_map) ;
        title('First estimation...', 'FontSize', FontSize) ;
        set(findall(gca,'-property','Interpreter'),'Interpreter','Latex') ;
        set(findall(gcf,'-property','Interpreter'),'Interpreter','Latex') ;
        plot_circle(fig_1, ...
            [list_y(l, :)', list_x(l, :)'], ...
            scale_plot*list_sig(l, :)', pix_cube, [1,1,1]) ;
        % Set marker size
        list_plot = get(gca(fig_1),'Children') ;
        for h = 1:length(list_plot)
            if strcmp('line', get(list_plot(h), 'Type'))
                set(list_plot(h), 'MarkerSize', MarkerSize, ...
                    'LineWidth', LineWidth) ;
            end
        end
        set(fig_1,'rend','painters','pos', [fig_pos, fig_size]) ;

        %% Simulating frame l
        for corner = 1:4
            % Extracting subpicture
            ind_y = pos2ind( [list_y(l, corner), list_x(l, corner)], ...
                pix_cube) ;
            ind_x = ind_y(2)+(-rad_ROI:rad_ROI) ;
            ind_y = ind_y(1)+(-rad_ROI:rad_ROI) ;
            y_frame = get_Fourier_vector(pix_cube.nb_y, 1) ;
            x_frame = get_Fourier_vector(pix_cube.nb_x, 1) ;
            y_corner = y_frame(ind_y) ;
            x_corner = x_frame(ind_x) ;

            % Simulating gaussian patterns
            Gauss = OpGauss(y_corner, x_corner, ...
                [], [], [], [], [], flag_norm) ;

            hyp_cube_sim(ind_y, ind_x, l) = Gauss * [ ...
                list_y(l, corner); list_x(l, corner); ...
                list_sig(l, corner); list_amp(l, corner); ...
                list_off(l)] ;
        end

        fig_2 = figure(2) ;
        imshow( hyp_cube_sim(:, :,l), color_scale) ;
        colormap(gca, color_map) ;
        title('First estimation...', 'FontSize', FontSize) ;
        set(findall(gca,'-property','Interpreter'),'Interpreter','Latex') ;
        set(findall(gcf,'-property','Interpreter'),'Interpreter','Latex') ;
        set(fig_2,'rend','painters','pos', [fig_pos, fig_size]) ;

        %% Asking the user if he is satisfied?
        str_user = input('Is this guess correct? (y/n)\n', 's') ;
        first_frame = ~strcmp(str_user, 'y') ;
    end

    %% Loop on the other wavelength
    for l = (l_first+1):min(l_end,pix_cube.nb_l)
        disp(['Analyzing wavelength ', num2str(l-l_first+1), '/', ...
            num2str((min(l_end,pix_cube.nb_l)-l_first+1))]) ;
        %% Initialization
        frame_l = hyp_cube(:,:,l) ;
        med_frame = median(frame_l(:)) ;
        list_off(l) = med_frame ;

        % Position of the previous wavelength
        list_y(l, :) = list_y(l-1, :) ;
        list_x(l, :) = list_x(l-1, :) ;
        for corner = 1:4
            ind_y = pos2ind( [list_y(l, corner), list_x(l, corner)], ...
                        pix_cube) ;
            ind_x = ind_y(2) ;
            ind_y = ind_y(1) ;
            list_amp(l, corner) = frame_l(ind_y, ind_x) ;
        end

        %% Refining the square
        % Loop on the corners
        if paralleling
            parfor corner = 1:4
                for cal = 1:nb_cal
                    % Refining position and amplitude
                    for cc = 1:length(list_cal_pat)
                        [list_y(l, corner), list_x(l, corner), ...
                            list_amp(l, corner)] = ...
                            find_pattern_pos_amp(frame_l, ...
                            [list_y(l, corner); list_x(l, corner); ...
                            list_amp(l, corner)], ...
                            rad_ROI, pix_cube, option_opti_amp_pos, ...
                            pattern_model, ...
                            list_sig(l, corner), med_frame, ...
                            list_cal_pat{cc}) ;
                    end

                    % Refining pattern parameter
                    list_sig(l, corner) = find_pattern_par(frame_l, ...
                        list_sig(l, corner), ...
                        rad_ROI, pix_cube, option_opti_par, ...
                        pattern_model, ...
                        list_y(l, corner), list_x(l, corner), ...
                        list_amp(l, corner), med_frame) ;
                end
            end
        else
            for corner = 1:4
                for cal = 1:nb_cal
                    % Refining position and amplitude
                    for cc = 1:length(list_cal_pat)
                        [list_y(l, corner), list_x(l, corner), ...
                            list_amp(l, corner)] = ...
                            find_pattern_pos_amp(frame_l, ...
                            [list_y(l, corner); list_x(l, corner); ...
                            list_amp(l, corner)], ...
                            rad_ROI, pix_cube, option_opti_amp_pos, ...
                            pattern_model, ...
                            list_sig(l, corner), med_frame, ...
                            list_cal_pat{cc}) ;
                    end

                    % Refining pattern parameter
                    list_sig(l, corner) = find_pattern_par(frame_l, ...
                        list_sig(l, corner), ...
                        rad_ROI, pix_cube, option_opti_par, ...
                        pattern_model, ...
                        list_y(l, corner), list_x(l, corner), ...
                        list_amp(l, corner), med_frame) ;
                end
            end
        end

        %% Simulating frame l
        for corner = 1:4
            % Extracting subpicture
            ind_y = pos2ind( [list_y(l, corner), list_x(l, corner)], ...
                pix_cube) ;
            ind_x = ind_y(2)+(-rad_ROI:rad_ROI) ;
            ind_y = ind_y(1)+(-rad_ROI:rad_ROI) ;
            y_frame = get_Fourier_vector(pix_cube.nb_y, 1) ;
            x_frame = get_Fourier_vector(pix_cube.nb_x, 1) ;
            y_corner = y_frame(ind_y) ;
            x_corner = x_frame(ind_x) ;

            % Simulating gaussian patterns
            Gauss = OpGauss(y_corner, x_corner, ...
                [], [], [], [], [], flag_norm) ;

            hyp_cube_sim(ind_y, ind_x, l) = Gauss * [ ...
                list_y(l, corner); list_x(l, corner); ...
                list_sig(l, corner); list_amp(l, corner); ...
                list_off(l)] ;
        end
    end

    %% Saving
    save_fits(hyp_cube_sim, 'hyp_cube_sim', save_path) ;
    save([save_path, 'list_y.mat'], 'list_y') ;
    save([save_path, 'list_x.mat'], 'list_x') ;
    save([save_path, 'list_sig.mat'], 'list_sig') ;
    save([save_path, 'list_amp.mat'], 'list_amp') ;
    save([save_path, 'list_off.mat'], 'list_off') ;
end

%% Estimating squares
if ~load_centroid
    % Loop on the wavelength
    for l = 1:pix_cube.nb_l
        disp(['Fitting wavelength ', num2str(l), '/', ...
            num2str(pix_cube.nb_l)]) ;
        %% Initialization
        % Centroid
        list_centroid(l,1) = mean(list_y(l,:)) ;
        list_centroid(l,2) = mean(list_x(l,:)) ;
        
        % Diagonal
        list_diag(l) = mean( ...
            ( ...
            (list_y(l,:)-list_centroid(l,1)).^2 + ...
            (list_x(l,:)-list_centroid(l,2)).^2).^0.5 ...
            ) ;
        
        % Orientation
        list_theta(l) = mean(mod( ...
            atand( ...
            (list_y(l,:)-list_centroid(l,1))./ ...
            (list_x(l,:)-list_centroid(l,2))), ...
            90)) ;
        
        %% Refinement of the square parameters
        dist2 = @(square_par)sum( ...
            ( ...
                square_par(1) + ...
                square_par(3)*sind(square_par(4)+(0:90:270)') - ...
                list_y(l,:)').^2 ...
            + ...
            ( ...
                square_par(2) + ...
                square_par(3)*cosd(square_par(4)+(0:90:270)') - ...
                list_x(l,:)').^2 ...
            ) ;
        square_par = run_OptiFMINSEARCH(dist2, [], [], option_opti, ...
            [list_centroid(l,:), list_diag(l), list_theta(l)]) ;
        
        %% Update the values
        % Centroid
        list_centroid(l,1) = square_par(1) ;
        list_centroid(l,2) = square_par(2) ;
        
        % Diagonal
        list_diag(l) = square_par(3) ;
        
        % Orientation
        list_theta(l) = square_par(4) ;
    end
    
    
    %% Saving
    save_fits(hyp_cube_sim, 'hyp_cube_sim', save_path) ;
    save([save_path, 'list_centroid.mat'], 'list_centroid') ;
    save([save_path, 'list_diag.mat'], 'list_diag') ;
    save([save_path, 'list_theta.mat'], 'list_theta') ;
end

%% Display results
center = [-Inf, -Inf] ;
flag_center = false(pix_cube.nb_l, 1) ;
flag_center(l_first:l_end) = true ;
center_new = median(list_centroid, 1) ;
while sum(center~=center_new)==2
    % Update center
    center = center_new ;
    
    % New estimation of the center
    center_new = mean(list_centroid(flag_center, :), 1) ;
    flag_center = sum((list_centroid-center_new).^2,2)<tol_center^2 ;
end

l = 1 ;
cont = true ;
fig_1 = figure(1) ;
while cont
    %% Initialization
    frame_l = hyp_cube(:,:,l) ;
    if flag_center(l) 
        color_l = 0.999*[1,1,1] ;
    else
        color_l = 0*[1,1,1] ;
    end
    
    % Display
    imshow(frame_l, color_scale) ;
    colormap(gca, color_map) ;
    title(['Fitting wavelength ', num2str(l),'...'], 'FontSize', FontSize) ;
    set(findall(gca,'-property','Interpreter'),'Interpreter','Latex') ;
    set(findall(gcf,'-property','Interpreter'),'Interpreter','Latex') ;
    plot_circle(fig_1, ...
        [list_y(l, :)', list_x(l, :)'], ...
        scale_plot*list_sig(l, :)', pix_cube, color_l) ;
    plot_rectangle(fig_1, ...
        [list_centroid(l, 1), list_centroid(l, 2)], ...
        2/sqrt(2)*list_diag(l)*[1,1], list_theta(l)-45, pix_cube, ...
        color_l) ;    
    
    % Set marker size
    list_plot = get(gca(fig_1),'Children') ;
    for h = 1:length(list_plot)
        if strcmp('line', get(list_plot(h), 'Type'))
            set(list_plot(h), 'MarkerSize', MarkerSize, ...
                'LineWidth', LineWidth) ;
        end
    end
    set(fig_1,'rend','painters','pos', [fig_pos, fig_size]) ;
    
    %% User interface
    waitforbuttonpress;
    value = double(get(gcf,'CurrentCharacter')) ;
    switch value
        case 27
            cont = false ;
        case {28, 30}
            l = l-1 ;
            l = max(1,l) ;
        case {29, 31}
            l = l+1 ;
            l = min(l,pix_cube.nb_l) ;
    end
end
save([save_path, 'pos_center.txt'], 'center', '-ascii', '-double') ;

%% Plot center
fig_1 = figure(1) ;
set(fig_1,'rend','painters','pos', [fig_pos, fig_size(2), fig_size(2)]) ;
clf(fig_1) ;
hold on
plot(list_centroid(~flag_center,1)-center(1), ...
    list_centroid(~flag_center,2)-center(2), '.red') ;
plot(list_centroid(flag_center,1)-center(1), ...
    list_centroid(flag_center,2)-center(2), '.green') ;
plot(0, 0, '.black') ;
hold off
% Set marker size
list_plot = get(gca(fig_1),'Children') ;
for h = 1:length(list_plot)
    if strcmp('line', get(list_plot(h), 'Type'))
        set(list_plot(h), 'MarkerSize', 1.5*MarkerSize, ...
            'LineWidth', LineWidth) ;
    end
end
axis equal ;
xticks(linspace(-nb_tol_plot*tol_center, nb_tol_plot*tol_center, 11)) ;
yticks(linspace(-nb_tol_plot*tol_center, nb_tol_plot*tol_center, 11)) ;
axis([-nb_tol_plot*tol_center, nb_tol_plot*tol_center, ...
    -nb_tol_plot*tol_center, nb_tol_plot*tol_center]) ;
set(gca,'FontSize',FontSize_axis);
title('Dispersion of the squares centers', 'FontSize', FontSize) ;
set(findall(gca,'-property','Interpreter'),'Interpreter','Latex') ;
set(findall(gcf,'-property','Interpreter'),'Interpreter','Latex') ;
if flag_png
    saveas(fig_1, [save_path, 'Dispersion of the squares centers'], ...
        'png') ;
end


%% Saving figures
if flag_png
    fig_1 = figure(1) ;
    for l = 1:pix_cube.nb_l
        %% Initialization
        frame_l = hyp_cube(:,:,l) ;
        if flag_center(l) 
            color_l = 0.999*[1,1,1] ;
        else
            color_l = 0*[1,1,1] ;
        end

        % Display
        imshow(frame_l, color_scale) ;
        colormap(gca, color_map) ;
        title(['Fitting wavelength ', num2str(l),'...'], 'FontSize', ...
            FontSize) ;
        set(findall(gca,'-property','Interpreter'),'Interpreter','Latex') ;
        set(findall(gcf,'-property','Interpreter'),'Interpreter','Latex') ;
        plot_circle(fig_1, ...
            [list_y(l, :)', list_x(l, :)'], ...
            scale_plot*list_sig(l, :)', pix_cube, color_l) ;
        plot_rectangle(fig_1, ...
            [list_centroid(l, 1), list_centroid(l, 2)], ...
            2/sqrt(2)*list_diag(l)*[1,1], list_theta(l)-45, pix_cube, ...
            color_l) ;    

        % Set marker size
        list_plot = get(gca(fig_1),'Children') ;
        for h = 1:length(list_plot)
            if strcmp('line', get(list_plot(h), 'Type'))
                set(list_plot(h), 'MarkerSize', MarkerSize, ...
                    'LineWidth', LineWidth) ;
            end
        end
        set(fig_1,'rend','painters','pos', [fig_pos, fig_size]) ;
        pause(0.1) ;
        saveas(gcf, [save_path, 'square_fit/', 'lambda_', ...
            num2str(l)], 'png') ;
    end
end

%% Analyze scaling
% Fitting coefficient
nb_flag = sum(flag_center) ;
list_lambda_flag = reshape(list_lambda(flag_center), [1, 1, nb_flag]) ;
list_diag_flag = reshape(list_diag(flag_center), [1, 1, nb_flag]) ;

[a, b] = fit_linear_law(list_lambda_flag, list_diag_flag, [], ...
    flag_offset) ;

fig_1 = figure(1) ;
clf(fig_1) ;
hold on
plot(list_lambda, a*list_lambda+b, '-b', 'LineWidth', LineWidth) ;
plot(list_lambda(flag_center), list_diag(flag_center), 'og', ...
    'LineWidth', LineWidth) ;
plot(list_lambda(~flag_center), list_diag(~flag_center), 'or', ...
    'LineWidth', LineWidth) ;
hold off
axis([list_lambda(1), list_lambda(pix_cube.nb_l), ...
    a*list_lambda(1)+b, a*list_lambda(pix_cube.nb_l)+b]) ;
title(['Wavelength scaling (a=',num2str(a), ' // b=', num2str(b), ')']) ;
set(fig_1,'rend','painters','pos', ...
    [fig_pos, fig_size(2), fig_size(2)/1.5]) ;
set(gca,'FontSize',FontSize_axis);
set(findall(gca,'-property','Interpreter'),'Interpreter','Latex') ;
xlabel('$\lambda$ ($\mu m$)', 'FontSize', FontSize) ;
ylabel('Diagonal (pixel)', 'FontSize', FontSize) ;
set(findall(gcf,'-property','Interpreter'),'Interpreter','Latex') ;

fig_2 = figure(2) ;
clf(fig_2) ;
hold on
plot(list_lambda, 0*list_lambda+b, '-b', 'LineWidth', LineWidth) ;
plot(list_lambda(flag_center), ...
    list_diag(flag_center)-a*list_lambda(flag_center), 'og', ...
    'LineWidth', LineWidth) ;
plot(list_lambda(~flag_center), ...
    list_diag(~flag_center)-a*list_lambda(~flag_center), ...
    'or', 'LineWidth', LineWidth) ;
hold off
axis([list_lambda(1), list_lambda(pix_cube.nb_l), ...
    b-nb_tol_plot*tol_center, b+nb_tol_plot*tol_center]) ;
title(['Wavelength scaling - residues (a=',num2str(a), ' // b=', ...
    num2str(b), ')']) ;
set(fig_2,'rend','painters','pos', ...
    [fig_pos+[500,0], fig_size(2), fig_size(2)/1.5]) ;
set(gca,'FontSize',FontSize_axis);
set(findall(gca,'-property','Interpreter'),'Interpreter','Latex') ;
xlabel('$\lambda$ ($\mu m$)', 'FontSize', FontSize) ;
ylabel('Diagonal (pixel)', 'FontSize', FontSize) ;
set(findall(gcf,'-property','Interpreter'),'Interpreter','Latex') ;

if flag_png
    saveas(fig_1, [save_path, 'Wavelength scaling'], 'png') ;
    saveas(fig_2, [save_path, 'Wavelength scaling - residues'], 'png') ;
end