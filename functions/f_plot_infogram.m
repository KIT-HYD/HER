function A = f_plot_infogram(her)
%% function to plot spatial characterization 

% -------------- Input --------------
% - her       struct  structure cointaining model definitions

% -------------- Version --------------
% - 2020/03/20 Stephanie Thiesen: intial version
% - 2020/03/23 Stephanie Thiesen: added infogram plot

% -------------- Script --------------
    dim_cal = length(her.mat_euc_distance_xy);
    
    % Delta_z_PMF x classes
    figure;                                                                                                                     
    plot(her.mat_euc_distance_xy, her.mat_diff_z, 'xb');                                                                                                                                    
    xlabel('Euclidean distance');
    ylabel('\Deltaz');
    hold on;
    %     xlim([0,hmax]);
    for i = 2 : length(her.edges_distance_classes)
        line([her.edges_distance_classes(i),her.edges_distance_classes(i)], get(gca, 'ylim'),'Color','red','LineStyle','--');
    end
    title({strcat('Infogram cloud / cal. set: ', num2str(dim_cal), ' obs.');''});
    pbaspect([1.5 1 1]);

    % Infogram
    hmax = max(her.mat_euc_distance_xy(:));
    figure;
    plot(her.bin_centers_distance_classes, her.H_diff_z_by_class,'Marker', '.', 'Color', 'red', 'MarkerSize', 20);
    title({strcat('Infogram / ', her.txt);''});
    pbaspect([1.5 1 1]);
    hold on;
    line([0 hmax],[her.H_diff_z her.H_diff_z], 'LineWidth', 2);
    for i = 2 : length(her.edges_distance_classes)
        line([her.edges_distance_classes(i),her.edges_distance_classes(i)], get(gca, 'ylim'),'Color','black','LineStyle','--');
    end
    legend({'Entropy of \Deltaz PMF by class', 'Entropy of the fullset \Deltaz PMF', 'Distance class limit'},'Location','southwest');
    xlabel('Euclidean distance');
    ylabel('Entropy [bit]');
    
    % # of pairs/class
    fig = figure();
    left_color = [0 0 0.4];
    right_color = [1 0 0];
    set(fig,'defaultAxesColorOrder',[left_color; right_color]);

    yyaxis right
    plot(1:her.n_lag, her.H_diff_z_by_class,'Marker', '.', 'Color', 'red', 'MarkerSize', 20);
    ylim([0, max(her.H_diff_z_by_class,[],'all')+max(her.H_diff_z_by_class,[],'all')/90])
    ylabel('Entropy [bit]');

    yyaxis left
    b = bar(1:her.n_lag,her.n_pairs_by_class);%, 'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5);
    b.FaceColor = 'flat';
    for i = 1:her.n_classes_range-1
        b.CData(i,:) = [0 0 .9];
    end
    xticks([1 5:5:her.n_lag]);
    ax = gca;
    ax.YMinorTick = 'on';
    ax.YAxis(1).MinorTickValues = 0:500:max(ylim);
    xlabel('Distance class');
    ylabel('Number of pairs');

    title({strcat('Infogram & Number of pairs (in \Deltaz histogram) by class');''});
    pbaspect([1.5 1 1]);
    str = {['# of observations: ' num2str(dim_cal) ' (cal. set)'], ['# of pairs: ' ...
        num2str(sum(her.n_pairs_by_class)) ' (' num2str(dim_cal) 'x' num2str(dim_cal-1) ')' ],...
        ['# of pairs into the range: ' num2str(sum(her.n_pairs_by_class(1:her.n_classes_range-1))) ...
        ' (' num2str(round(100*sum(her.n_pairs_by_class(1:her.n_classes_range-1))/sum(her.n_pairs_by_class),1)) '%)']};% ['# of pairs out of the range: ' num2str(sum(her.n_pairs_by_class(her.n_classes_range:end))) ' (' num2str(round(100*sum(her.n_pairs_by_class(her.n_classes_range:end))/sum(her.n_pairs_by_class),1)) '%)']};
    text(her.n_lag/1.3,max(her.n_pairs_by_class)/1,str);

    % delta-z PMF/class
    figure;
    ncols = 5;
    nrows = ceil(((length(her.bin_centers_distance_classes_range)) / ncols));
    for i = 1:her.n_classes_range
        subplot(nrows,ncols,i);
        bar(her.bin_centers_edges_diff_z, her.pmf_diff_z_by_class_obs_range(i,:));%, 0.9, 'histc');
        title(['class ' num2str(i) ' / dist.: (' num2str(round(her.edges_distance_classes_range(i),2)) ', ' num2str(floor(her.edges_distance_classes_range(i+1))) '] ']);
        pbaspect([3 1 1]);
        xlabel('\Deltaz');
        ylim([0, max(her.pmf_diff_z_by_class_obs_range,[],'all')+max(her.pmf_diff_z_by_class_obs_range,[],'all')/50]);
        xlim([-5 5]);
    end
    title(['class ' num2str(i) ' (full dataset PMF)']);
    sgtitle({strcat('\Deltaz PMF by class');''});

    % EXTRA step: Variogram
    lag_variance = NaN(length(her.edges_distance_classes)-1,1); %vector for the semivariance of the observations
    n_obs_by_lag = zeros(length(her.edges_distance_classes)-1,1); %vector for counting the number of pair points in each lag class
    mat_square_diff_z = her.mat_diff_z.^2; %matrix of the squared semivariance
    for i = 1 : (length(her.edges_distance_classes)-1) %for each distance class
        idx = her.mat_euc_distance_xy > her.edges_distance_classes(i) & her.mat_euc_distance_xy <= her.edges_distance_classes(i+1); %find the observations within the current lag class
        n_obs_by_lag(i)= sum(idx(:)); %count the number of values
        lag_variance(i) = 1/2 * mean(mat_square_diff_z(idx)); %calculate the semivariance
    end

    %Infogram + Variogram
    fig = figure();
    left_color = [0 0 0.4];
    right_color = [.7 0 0];
    set(fig,'defaultAxesColorOrder',[left_color; right_color]);

    plot(her.bin_centers_distance_classes, her.H_diff_z_by_class./max(her.H_diff_z_by_class),'Marker', '.', 'MarkerSize', 20);

    yyaxis left
    ylabel('Normalized entropy [bit]');
    xlim([0 her.bin_centers_distance_classes(end-1)]);
    ylim([0 1]);

    yyaxis right
    plot(her.bin_centers_distance_classes, lag_variance./(max(lag_variance)),'Marker', '.', 'MarkerSize', 20);
    xlabel('lag');
    ylabel('Normalized dissimilarity (1/2*mean(\Deltaz²))');
    hold on 
    for i = 2 : length(her.edges_distance_classes)
        line([her.edges_distance_classes(i),her.edges_distance_classes(i)], get(gca, 'ylim'),'Color','black','LineStyle','--');
    end
    ylim([0 1]);
    title({'Normalized Infogram and Experimental variogram';''});


end