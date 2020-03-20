function [z_target_entropy_pred_plot, z_target_mean_pred_GRID_plot, z_target_median_pred_GRID_plot, z_target_mode_pred_GRID_plot, varargout] = f_extract_pmf_statistics(x_target_grid, y_target_grid, pmf_pred_nn, bin_centers_edges_z, varargin)
%% function to extract PMF entropy, mean, median, mode, exceed probability

% -------------- Input --------------
% - x_target_grid       [T,1]   x coordidates of the target set
% - y_target_grid       [T,1]   x coordidates of the target set
% - pmf_pred_nn         {1,T}   predicted z PMF for targets 
% - bin_centers_edges_z [1,n]   bin centers of the z PMF
% - varargin            h       threshold of z

% -------------- Output --------------
% - z_target_entropy_pred_plot     [T,T]   entropy of the predicted target z_PMFs
% - z_target_mean_pred_GRID_plot   [T,T]   mean of the predicted target z_PMFs
% - z_target_median_pred_GRID_plot [T,T]   median of the predicted target z_PMFs
% - z_target_mode_pred_GRID_plot   [T,T]   mode of the predicted target z_PMFs
% - varargout                      [T,T]   probability of predicted target z_PMFs
%                                          exceeding varargin (threshold of z)                                        of the 

% -------------- Version --------------
% - 2020/03/20 Stephanie Thiesen: intial version

% -------------- Script --------------
    %calculate the entropy of the z PMFs for all predicted points 
    H_z_pmf_by_class_pred = NaN(numel(pmf_pred_nn),1); %vector for the entropy of the grid (xi,yi_analized)
    for target_ = 1:numel(pmf_pred_nn) %for each target
        H_z_pmf_by_class_pred(target_,1) = f_entropy(cell2mat(pmf_pred_nn(1,target_))); %calculate the entropy
    end
    z_target_entropy_pred_plot = reshape(H_z_pmf_by_class_pred, size(y_target_grid,1), size(y_target_grid,2));

    % z mean, median, mode,
    z_target_mean_pred = NaN(1,numel(pmf_pred_nn));
    z_target_median_pred = NaN(1,numel(pmf_pred_nn));
    z_target_mode_pred = NaN(1,numel(pmf_pred_nn));
    z_target_probability_pred = NaN(1,numel(pmf_pred_nn));
    
    for target_ = 1:numel(pmf_pred_nn) %xi,yi_analized 
        %mean
        z_target_mean_pred(1,target_) = (sum(bin_centers_edges_z .* (pmf_pred_nn{1,target_}))) / sum(pmf_pred_nn{1,target_}); %interpolated mean (not the bin center)
        %median
        cmf_median_ = cumsum(pmf_pred_nn{1,target_});
        idx_median_ = find(cmf_median_ >= 0.5, 1);
        x_median_ = [bin_centers_edges_z(idx_median_-1) bin_centers_edges_z(idx_median_) bin_centers_edges_z(idx_median_+1)];
        y_median_ = [cmf_median_(idx_median_-1) cmf_median_(idx_median_) cmf_median_(idx_median_+1)];
        z_target_median_pred(1,target_) = interp1(y_median_,x_median_,0.5); %interpolated median (not the bin center)
        %mode
        [~,idx_] = max(cell2mat(pmf_pred_nn(1,target_)));
        z_target_mode_pred(1,target_) = bin_centers_edges_z(idx_); %bin center
    end
    if ~isnan(varargin{1})
        if length(varargin) >= 1 %nargin >= 4
            thres = varargin{1};
            for target_ = 1:numel(pmf_pred_nn) %xi,yi_analized 
                %probability of exceeding thres
                bin_thres_ = sum(bin_centers_edges_z <= thres);
                if bin_centers_edges_z(bin_thres_) < thres
                    bin_thres_ = bin_thres_+1;
                end 
                cmf_probab = cumsum(pmf_pred_nn{1,target_});
                z_target_probability_pred(1,target_) = 1 - cmf_probab(bin_thres_);
            end
        end
            z_target_probability_pred_GRID_plot = reshape(z_target_probability_pred, size(y_target_grid,1), size(y_target_grid,2));
    else
        z_target_probability_pred_GRID_plot = NaN;
    end

    z_target_mean_pred_GRID_plot = reshape(z_target_mean_pred, size(y_target_grid,1), size(y_target_grid,2));
    z_target_median_pred_GRID_plot = reshape(z_target_median_pred, size(y_target_grid,1), size(y_target_grid,2));
    z_target_mode_pred_GRID_plot = reshape(z_target_mode_pred, size(y_target_grid,1), size(y_target_grid,2));

    if nargout >= 2
        varargout{1} = z_target_probability_pred_GRID_plot;
    end

end