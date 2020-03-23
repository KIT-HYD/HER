function A = f_plot_prediction(z_plot, z_entropy_plot, pmf_pred, x_target, y_target, x_target_grid, y_target_grid, x, y, z, idx_cal, her, x_pmf, y_pmf, true_on, varargin)
%% function to plot HER prediction and original dataset
% -------------- Input -------------- 
% - z_plot              [T,T]          z results to be ploted
% - z_entropy_plot      [T,T]          entropy results to be ploted
% - pmf_pred            {1,T}          z_PMFs of the predicted locations
% - x_target_grid; y_target_grid   [T,T]    x,y coordinates of the GRID
% - x; y                [n,1]          x,y coordinates of the original dataset
% - z                   [n,1]          z true values the predicted locations
% - idx_cal             [1,c]          index of the calibration set
% - x_pmf; y_pmf        [1,p]          x,y coordinates to plot predicted PMF
% - true_on             char   'N' if no true z is know for x_pmf and y_pmf  
%                              'Y'  if there is a true z for x_pmf and y_pmf 
% - varargin (shp_basin)  struc         basin shapefile

% -------------- Version --------------
% - 2020/03/20 Stephanie Thiesen: intial version
% - 2020/03/23 Stephanie Thiesen: removed possibility to plot without GRID

% -------------- Script --------------
    if length(varargin) >= 1
        shp_basin = varargin{1};
    end

    clims = [min([z(:);z_plot(:)]) max([z(:);z_plot(:)])];
    sz = 10;
    for i = 1:length(x_pmf)
        randplot_set(i,1) = find(and(ismember(x_target, x_pmf(i)),ismember(y_target, y_pmf(i))));
    end   
    ncols = 2;
    nrows = ceil((length(x_pmf) / ncols));

        figure;
        ax1 = subplot(1,3,1);
        pcolor(x_target_grid,y_target_grid,z);
        hold on
        scatter(x(idx_cal), y(idx_cal), 1+70*normalize(z(idx_cal),'range'),'k+','LineWidth',1);
        scatter(x_pmf, y_pmf,100,'s','MarkerEdgeColor', [0.9 0.9 0.9], 'LineWidth',1.5);
        colormap(gca,parula(6));
        title({'Original';''});
        xlabel('x');
        ylabel('y');
        h=colorbar;
        caxis(clims);
        shading flat;
        pbaspect([1 1 1]);
        if exist('shp_basin','var')
            symspec_ = makesymbolspec('Polygon', ...
               {'ROCK', 'basin','FaceColor', [1 1 1], 'FaceAlpha',0,'LineWidth',2, 'EdgeColor', [0.9 0.9 0.9]});
            mapshow(shp_basin,'SymbolSpec', symspec_);
        end

        ax2 = subplot(1,3,2);
        hold on
        pcolor(x_target_grid, y_target_grid, z_plot);
        scatter(x(idx_cal), y(idx_cal), 1+70*normalize(z(idx_cal),'range'),'k+','LineWidth',1);
        scatter(x_pmf, y_pmf,100,'s','MarkerEdgeColor', [0.9 0.9 0.9], 'LineWidth',1.5);
        colormap(gca,parula(6));
        title({'HER E-type';''});
        h=colorbar;
        caxis(clims);
        shading flat;
        pbaspect([1 1 1]);
        if exist('shp_basin','var')
            mapshow(shp_basin,'SymbolSpec', symspec_);
        end

        % plot entropy map
        ax3 = subplot(1,3,3);
        clims_h = [min(z_entropy_plot(:)) max(z_entropy_plot(:))];
        hold on
        pcolor(x_target_grid, y_target_grid, z_entropy_plot);
        scatter(x(idx_cal), y(idx_cal), 1+70*normalize(z(idx_cal),'range'),'k+','LineWidth',1);
        scatter(x_pmf, y_pmf,100,'s','MarkerEdgeColor', [0.9 0.9 0.9], 'LineWidth',1.5);
        colormap(gca,spring(8));
        title({'Entropy map [bits]';''});
        shading flat;
        % set(gca,'ColorScale','log')
        h=colorbar;
        hold on;
        scatter(x(idx_cal), y(idx_cal), 1+70*normalize(z(idx_cal),'range'),'k+','LineWidth',1);
        caxis(clims_h);
        pbaspect([1 1 1]);
        if exist('shp_basin','var')
            symspec_ = makesymbolspec('Polygon', ...
               {'ROCK', 'basin','FaceColor', [1 1 1], 'FaceAlpha',0,'LineWidth',2, 'EdgeColor', [0.9 0.9 0.9]});
            mapshow(shp_basin,'SymbolSpec', symspec_);
        end
        sgtitle({strcat('Maps | Dataset: ', her.txt);''});
        linkaxes([ax1,ax2, ax3],'xy')
        legend(['H'],['cal'],['selected PMF'],'Location','northwest');

        %PMF
        figure;
        i = 1;
        if true_on == 'N'
            for target = randplot_set'
                subplot(nrows,ncols,i);
                h1 = bar(her.bin_centers_edges_z, cell2mat(pmf_pred(1,target)), 'k');
                xlim([-10 10]);
                ylim([0 0.2]);
                xlabel('z');
                title({['target ' num2str(target)] ,['x=' num2str(x_target(target)) '; y=' num2str(y_target(target)) ] ,['H = ' num2str(round(f_entropy(cell2mat(pmf_pred(1,target))),2)) ' bits']});
                i=i+1;
                pbaspect([1 1 1]);
            end
        end
        if true_on == 'Y'
             for target = randplot_set'
                subplot(nrows,ncols,i);
                h1 = bar(her.bin_centers_edges_z, cell2mat(pmf_pred(1,target)), 'k');
                hold on;
                h2 = plot([z(target) z(target)], [0, 0.5], 'm--', 'LineWidth',1.5);
                xlim([-10 10]);
                ylim([0 0.2]);
                xlabel('z');
                title({['target ' num2str(target)] ,['x=' num2str(x_target(target)) '; y=' num2str(y_target(target)) ] ,['H = ' num2str(round(f_entropy(cell2mat(pmf_pred(1,target))),2)) ' bits']});
                i=i+1;
                pbaspect([1 1 1]);
             end            
        end
        set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        if true_on == 'Y'
            legend({['true z']});
        end
        sgtitle({strcat('Predicted z PMF / ', her.txt);''});


end