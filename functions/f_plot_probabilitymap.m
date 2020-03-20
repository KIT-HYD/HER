function A = f_plot_probabilitymap(z_prob_plot, z_thresh, txt, x_target, y_target, x_target_grid, y_target_grid, x, y, z, idx_cal, idx_val, varargin)
%% function to plot the probabilty of z > z_thresh (HER probability map)
% -------------- Input -------------- 
% - z_prob_plot         [T,T]or[T,1]   probability results to be ploted
% - z_thresh             t             z threshold 
% - txt                 char           dataset name
% - x_target; y_target  [T,T]or[T,1]   x,y coordinates of the predicted values
% - x_target_grid; y_target_grid       GRID based on the original dataset
% - x; y                [n,1]          x,y coordinates of the original dataset
% - z                   [n,1]          z true values the predicted locations
% - idx_cal             [1,c]          index of the calibration set
% - idx_val             [1,v]          index of the validation set
% - varargin (shp_basin)  struc        basin shapefile

% -------------- Version --------------
% - 2020/03/20 Stephanie Thiesen: intial version

% -------------- Script --------------

    if length(varargin) >= 1
        shp_basin = varargin{1};
    end
    
    vq_pred = griddata(x_target,y_target,z_prob_plot(:),x_target_grid,y_target_grid);
    % probabilty map
    figure;
    pcolor(x_target_grid, y_target_grid, vq_pred);

    xlabel('x');
    ylabel('y');
    shading flat;
    colormap(flipud(gray(60)));
    % set(gca)
    h=colorbar;
    hold on;
    scatter(x(idx_cal), y(idx_cal), 1+70*normalize(z(idx_cal),'range'),'r+','LineWidth',1);
    scatter(x(idx_val), y(idx_val),1+70*normalize(z(idx_val),'range'),'g+','LineWidth',1);
    caxis([0 1]);
    title({strcat('Probability map [Probability Z>',num2str(z_thresh),'] / ', txt);''});
    pbaspect([1 1 1]);
    limits = [0 15 30 45 60 75 100];
    colormap([repmat([.9 .9 .9],limits(2)-limits(1)+1,1);repmat([0.75 0.75 0.75],limits(3)-limits(2)+1,1); repmat([0.6 0.6 0.6],limits(4)-limits(3)+1,1);repmat([0.4 0.4 0.4],limits(5)-limits(4)+1,1); repmat([0.3 0.3 0.3],limits(6)-limits(5)+1,1); repmat([0.1 0.1 0.1],limits(7)-limits(6)+1,1)]);
    if exist('shp_basin','var')
        symspec_ = makesymbolspec('Polygon', ...
           {'ROCK', 'basin','FaceColor', [1 1 1], 'FaceAlpha',0,'LineWidth',2, 'EdgeColor', [0.2 0.4 1]});
        hold on
        mapshow(shp_basin,'SymbolSpec', symspec_);
    end
    legend(strcat('Probability Z>',num2str(z_thresh)), 'cal','val');

end