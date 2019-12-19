%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to plot the fitted law
% 
% Created: 03/04/2019 (mm/dd/yyyy)
% Anthony Berdeu (Laboratoire Hubert Curien)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Input
%%%%%%%%
% #list_exp_time# list of the exposure times
%
% #list_adu# list of the value of adu
%
% #list_fit# list of the fitted values (in the order of the exposure times)
%
% #list_var# list of the estimated variances
%
% #list_var_MAD# list of the robust estimation of the variance on the
% variance
%
% #flag_fig# number of the first figure
%
% #fig_pos# position of the figure
%
% #fig_size# size of the figures
%
% #LineWidth# Width of the plotted lines
%
% #MarkerSize# Size of the markers
%
% #[FontSize_axis, FontSize]# font size of the axes and legends
%
% #color# color of the plot
%
%%%%%%%%
% Ouput
%%%%%%%%
% #[fig_1, fig_2]# Handles to the figures
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig_1, fig_2] = ...
    plot_linear_law(list_exp_time, list_adu, list_fit, list_var, ...
    list_var_MAD, flag_fig, fig_pos, fig_size, LineWidth, MarkerSize, ...
    FontSize_axis, FontSize, color)

    %% Linear law with the median dark flux
    % New order
    [~,list_ind] = sort(list_adu) ;
    
    fig_1 = figure(flag_fig) ;
    clf(fig_1) ;
    set(fig_1,'rend','painters','pos', [fig_pos, fig_size]) ;
    hold on
    errorbar(list_adu, list_var, list_var_MAD, ...
        'o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize, ...
        'Color', color) ;
    plot(list_adu(list_ind), ...
        list_fit(list_ind), ...
        'LineWidth', LineWidth, 'MarkerSize', MarkerSize, ...
        'Color', color) ;
    hold off
    axis([0, 10*floor(max(list_adu)/10+1), ...
        0, max(list_var+list_var_MAD)]) ;
    set(gca,'FontSize',FontSize_axis);
    set(findall(gca,'-property','Interpreter'),'Interpreter','Latex') ;
    xlabel('Median flux (adu)', 'FontSize', FontSize) ;
    ylabel('Variance (adu)', 'FontSize', FontSize) ;
    set(findall(gcf,'-property','Interpreter'),'Interpreter','Latex') ;

    %% Linear law with the exposure time
    fig_2 = figure(flag_fig+1) ;
    clf(fig_2) ;
    set(fig_2,'rend','painters','pos', [fig_pos, fig_size]) ;
    hold on
    errorbar(list_exp_time, list_var,list_var_MAD, ...
        'o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize, ...
        'Color', color) ;
    plot(list_exp_time, list_fit, ...
        'LineWidth', LineWidth, 'MarkerSize', MarkerSize, ...
        'Color', color) ;
    hold off
    axis([0, 10*floor(max(list_exp_time)/10+1), ...
        0, max(list_var+list_var_MAD)]) ;
    set(gca,'FontSize',FontSize_axis);
    set(findall(gca,'-property','Interpreter'),'Interpreter','Latex') ;
    xlabel('Exposure time (s)', 'FontSize', FontSize) ;
    ylabel('Variance (adu)', 'FontSize', FontSize) ;
    set(findall(gcf,'-property','Interpreter'),'Interpreter','Latex') ;
end