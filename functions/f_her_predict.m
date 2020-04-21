function [pmf_pred_nn, target_idx_zero_neigh_pred] = f_her_predict(x_cal, y_cal,z_cal, x_target, y_target, her)
%% function to predict z PMFs using HER

% -------------- Input --------------
% - x_cal     [N,1]   x coordidates of the calibration set
% - y_cal     [N,1]   y coordidates of the calibration set
% - z_cal     [N,1]   z values of the calibration set (variable under study)
% - x_target  [T,1]   x coordidates of the target set
% - y_target  [T,1]   y coordidates of the target set
% - her       struct  structure cointaining model definitions
%   her.aggregation_type  char  definition of the aggregation method     
%                         'and'   log linear aggregation (f_loglinear_aggregation)
%                         'or'    linear aggregation (f_linear_aggregation)
%                         'andor' log linear aggregation between 'and' & 'or'
%       her.n_neighbor        nn    maximum number of neighbord to be 
%                                   considered for interpolation
%       her.n_classes_range   r     number of distance classes inside the 
%                                   range + 1 (full dataset class)
%       her.edges_z          [1,e]  edges of the bins of the z_PMF 
%       her.edges_distance_classes_range  [1,r+1]  edges of distance classes 
%                                                  until range + max. distance class
%   - variables automatically defined -
%       her.n_bins_z
%       her.edges_diff_z_shift
%       her.n_bins_shift
%       her.n_bins_z_per_bin_shift
%       her.binwidth_shift
%       her.pmf_diff_z_by_class_obs_range_shift
 
% -------------- Output --------------
% - pmf_pred_nn                 {1,T}   predicted z PMF for targets
% - target_idx_zero_neigh_pred  [1,i]   extract target with neighbors only beyond the range

% -------------- Version --------------
% - 2020/03/20 Stephanie Thiesen: intial version

% -------------- Script --------------
% Step 1: Distance of the target to its neighbors
    % Calculate the euclidean distance between target (column) and their 
    % known neighbors/observation (row)
    mat_euc_distance_obs_target_pred = NaN(length(x_cal),length(x_target(:))); % create empty vector
    for target = 1:numel(x_target)
        for obs = 1:length(x_cal)
            mat_euc_distance_obs_target_pred(obs,target) = f_euclidean_dist(x_target(target), y_target(target), x_cal(obs), y_cal(obs));
        end
    end

% Step 2: Identify the class of each neighbor
    % Classify the observations according to their distance to the target  
    classes_obs_target_pred = zeros(length(x_cal),length(x_target));
    for target = 1:numel(x_target)
        for obs = 1:length(x_cal)
            for class_ = 1 : her.n_classes_range
                if mat_euc_distance_obs_target_pred(obs,target) > her.edges_distance_classes_range(class_) & mat_euc_distance_obs_target_pred(obs,target) <= her.edges_distance_classes_range(class_+1); % in case it is within the current lag class
                    classes_obs_target_pred(obs,target) = class_; 
                end          
            end 
            if mat_euc_distance_obs_target_pred(obs,target) == 0 %if the distance between points is zero, it will receive the first class identifier
               classes_obs_target_pred(obs,target) = 1; 
            end
            if mat_euc_distance_obs_target_pred(obs,target) > her.edges_distance_classes_range(class_+1)               
                classes_obs_target_pred(obs,target) = her.n_classes_range(end); 
            end
        end    
    end  
    % Identification of the neighborless observations
    [n_obs_by_class_pred, edgge_class] = hist(classes_obs_target_pred(:,:), -0.5:1:her.n_classes_range); %getting the # of observations by class for each target (row)
    target_idx_zero_neigh_pred = [];
    for target = 1 : length(x_target) %for each target
        if n_obs_by_class_pred(end,target) == sum(n_obs_by_class_pred(2:end,target))
            target_idx_zero_neigh_pred = [target_idx_zero_neigh_pred target];
        end
    end
    % Finding the closest neighbors
    idx_nn = NaN(her.n_neighbor,length(z_cal));
    mat_euc_distance_obs_target_pred_nn = NaN(her.n_neighbor,length(z_cal));
    classes_obs_target_pred_nn = NaN(her.n_neighbor,length(z_cal));
    z_cal_nn = NaN(her.n_neighbor,length(z_cal));
    for target = 1 : length(x_target) %for each target identidy the closest neighbors
        [~,idx_] = sort(mat_euc_distance_obs_target_pred(:,target));
        idx_nn(:, target) = idx_(1:her.n_neighbor);
        mat_euc_distance_obs_target_pred_nn(:,target) = mat_euc_distance_obs_target_pred(idx_nn(:, target),target);
        classes_obs_target_pred_nn(:,target) = classes_obs_target_pred(idx_nn(:, target),target);
        z_cal_nn(:,target) = z_cal(idx_nn(:, target));
    end

% Step 3: Normalized weights
    weights_or_pred_nn = zeros(her.n_neighbor,length(x_target));
    normalized_weight_or_pred_nn = zeros(her.n_neighbor,length(x_target));
    weights_or_cont_pred_nn = zeros(her.n_neighbor,length(x_target));
    normalized_weights_or_cont_pred_nn = zeros(her.n_neighbor,length(x_target));
    weights_and_pred_nn = zeros(her.n_neighbor,length(x_target));
    weights_and_cont_pred_nn = zeros(her.n_neighbor,length(x_target));
    for target = 1:length(x_target)
        for obs = 1:her.n_neighbor
            weights_or_pred_nn(obs,target) = her.best_w_OR(classes_obs_target_pred_nn(obs,target)); %discrete OR weights
            weights_or_cont_pred_nn(obs,target) =  her.best_w_OR(classes_obs_target_pred_nn(obs,target)) + her.w_OR_slope(classes_obs_target_pred_nn(obs,target)) * ( mat_euc_distance_obs_target_pred_nn(obs,target) - her.edges_distance_classes_range(classes_obs_target_pred_nn(obs,target)) ); %continuous OR weights
            weights_and_pred_nn(obs,target) = her.best_w_AND(classes_obs_target_pred_nn(obs,target)); %discrete OR weights
            weights_and_cont_pred_nn(obs,target) =  her.best_w_AND(classes_obs_target_pred_nn(obs,target)) + her.w_AND_slope(classes_obs_target_pred_nn(obs,target)) * ( mat_euc_distance_obs_target_pred_nn(obs,target) - her.edges_distance_classes_range(classes_obs_target_pred_nn(obs,target)) ); %continuous OR weights
        end
    end

    for target = 1:length(x_target)
        for obs = 1:her.n_neighbor
            normalized_weight_or_pred_nn(obs,target) = weights_or_pred_nn(obs,target) ./ sum(weights_or_pred_nn(:,target)); %discrete OR weights
            normalized_weights_or_cont_pred_nn(obs,target) = weights_or_cont_pred_nn(obs,target) ./ sum(weights_or_cont_pred_nn(:,target)); %continuous OR weights
         end
    end

% Steps 4: Obtaining the z PMF of the known observations(obs_diff_z + z) & z PMF with
% aggregation
    pmf_pred_nn = cell(1,numel(x_target));
    target_idx_without_AND_pred = []; %vector for the observations which do not support AND combination
    i= 0;
    pmf_z_target_given_neigh_pred_ = nan(her.n_neighbor,her.n_bins_z); %z PMF for each neighbor (row)
    diff_n_bins_shift = her.n_bins_z * her.n_bins_z_per_bin_shift - her.n_bins_shift; %how many small bins we will need to add (before and/or after) to have "numbins_z_per_bin_shift" bins for each bin in the final z PMF
    for target = 1:numel(x_target)
        % z PMF of the known observations(obs_diff_z + z)
        for obs = 1:her.n_neighbor
            class_ = classes_obs_target_pred_nn(obs,target); 
            bins_shift = round((her.edges_diff_z_shift(1) + z_cal_nn(obs,target) - her.edges_z(1)) / her.binwidth_shift);
            idx_bins_shift = reshape( ...           % index matrix to find which bins from the shift PMF go into each bin of the z PMF
                  [ones(1,bins_shift) 2:her.n_bins_shift+1 ones(1,int16(diff_n_bins_shift) - bins_shift)], ... %fill edges with 1 (zero bin)
                  int16(her.n_bins_z_per_bin_shift), her.n_bins_z); % matrix with one column per z PMF bin and "numbins_z_per_bin_shift" rows
            pmf_ = her.pmf_diff_z_by_class_obs_range_shift(class_, :);
            pmf_z_target_given_neigh_pred_(obs,:) = sum(pmf_(idx_bins_shift),1); %group small bins from "her.pmf_diff_z_by_class_obs_range_shift" into large z PMF bins
            pmf_z_target_given_neigh_pred_(obs,:) = pmf_z_target_given_neigh_pred_(obs,:) / sum(pmf_z_target_given_neigh_pred_(obs,:)); %normalize
        end
        % predicting z PMF with aggregation 
        idx = [1:her.n_neighbor];
        pmfs_ = pmf_z_target_given_neigh_pred_(idx,:);
        if strcmp(her.aggregation_type,'andor')
            weights_or_ = normalized_weights_or_cont_pred_nn(idx,target);
            weights_and_ = weights_and_cont_pred_nn(idx,target);
            pmf_and_ = f_loglinear_aggregation(pmfs_,weights_and_);
            pmf_or_ = f_linear_aggregation(pmfs_,weights_or_);
            pmf_ = f_loglinear_aggregation([pmf_and_;pmf_or_],[her.best_alpha;her.best_beta]);
            [ pmf_pred_nn{1,target} ] = pmf_;
        elseif strcmp(her.aggregation_type,'and')
            weights_and_ = weights_and_cont_pred_nn(idx,target);
            pmf_and_ = f_loglinear_aggregation(pmfs_,weights_and_);
            [ pmf_pred_nn{1,target} ] = pmf_and_;
        elseif strcmp(her.aggregation_type,'or')
            weights_or_ = normalized_weights_or_cont_pred_nn(idx,target);
            pmf_or_ = f_linear_aggregation(pmfs_,weights_or_);
            [ pmf_pred_nn{1,target} ] = pmf_or_;
        end   
    end
end
