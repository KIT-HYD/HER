function A = f_plot_weights(her)
%% function to plot optimum weights

% -------------- Input --------------
% - her       struct  structure cointaining model definitions

% -------------- Version --------------
% - 2020/03/20 Stephanie Thiesen: intial version

% -------------- Script --------------
    figure;
    hold on;
    plot(her.edges_distance_classes_range(1:length(her.edges_distance_classes_range)-1), her.best_w_OR, 'Marker', '.','MarkerSize', 15, 'LineWidth', 2 );
    plot(her.edges_distance_classes_range(1:length(her.edges_distance_classes_range)-1), her.best_w_AND, 'Marker', '.','MarkerSize', 15, 'LineWidth', 2 );
    ylim([0 1]);

    for i = 2 : length(her.edges_distance_classes_range)
        line([her.edges_distance_classes_range(i-1),her.edges_distance_classes_range(i-1)], get(gca, 'ylim'),'Color','black','LineStyle','--');
    end
    title({strcat('Optimum weights by class');''});
    pbaspect([1.5 1 1]);
    xlabel('Euclidean distance');
    ylabel('Optimum weights');
    legend({[strcat('OR weights (DKL_O_R: ', num2str(round(her.DKL_w_OR,2)),')')]...
        [strcat('AND weights (DKL_A_N_D: ', num2str(round(her.DKL_w_AND,2)),')')]});

    figure;
    X = categorical({'Beta','Alpha'});
    X = reordercats(X,{'Beta','Alpha'});
    b = bar(X,[her.best_beta her.best_alpha]);
    b.FaceColor = 'flat';
    b.CData(2,:) = [.5 0 .5];
    ylabel('Optimum weights');
    ylim([0 1.2]);
    pbaspect([1.5 1 1]);
    title({strcat('Optimum weights (log-linear combination of OR and AND) | DKL_\alpha_\beta: ', num2str(round(her.DKL_w_alpha_beta,2)));''});
end   