function [her] = f_her_weight(x_cal, y_cal, z_cal, her)
%% function to optimize the weights of the classes

% -------------- Input --------------
% - x_cal     [N,1]   x coordidates of the calibration set
% - y_cal     [N,1]   y coordidates of the calibration set
% - z_cal     [N,1]   z values of the calibration set (variable under study)
% - her       struct  structure cointaining model definitions

% -------------- Output --------------
% - her       struct  structure cointaining model definitions

% -------------- Version --------------
% - 2020/04/10 Stephanie Thiesen: intial version

% -------------- Script --------------

% Geo2: Weight optimization 
% (leave-one-out cross validation on calibration dataset)

% step 2: calculate the distance of the target to its neighbors 
% calculate the euclidean distance between target and its neighbors(observations)
    % target: point to be predicted (left-out observation to be predicted and validated)
    % neighbors or observation: observation from the calibration set
mat_euc_distance_obs_target_opt = NaN(length(z_cal),length(z_cal)); %matrix for the distances (row:neighbors; column: target) 
for target = 1 : length(z_cal) %for each target
    for i = 1:length(z_cal) %for each neighbor
        mat_euc_distance_obs_target_opt(i,target) = f_euclidean_dist(x_cal(target), y_cal(target), x_cal(i), y_cal(i)); %calculate the distance between the neighbor and the target
    end
end 

% step 3: classify the neighbors according to the target
%according to its distance to the target, indentify the class of each neighbor 
classes_obs_target_opt = zeros(length(z_cal),length(z_cal)); %matrix for the neighbor class (row:neighbors; column: target). Class zero when the neighbor = target
for target = 1 : length(z_cal) %for each target
    for i = 1:length(z_cal) %for each neighbor
        for class_ = 1 : her.n_classes_range %for each class
            if mat_euc_distance_obs_target_opt(i,target) > her.edges_distance_classes_range(class_) & mat_euc_distance_obs_target_opt(i,target) <= her.edges_distance_classes_range(class_+1) % in case it is within the current lag class
                classes_obs_target_opt(i,target) = class_; %save the identified class
            end          
        end
    end    
end

% step 4: neighbor cardinality calculation
% plot the mean of the # of measurement neighbors by class
her.nmean_pairs_by_class = NaN(her.n_classes_range,1);
[her.n_pairs_by_class_obs, her.edges_n_obs] = hist(classes_obs_target_opt(:,:), -0.5:1:her.n_classes_range); %getting the # of neighbors by class (row) for each target (column)
for class_ = 1 : her.n_classes_range+1 %for each class (into the range + out of the range)
    her.nmean_pairs_by_class(class_) = mean(her.n_pairs_by_class_obs(class_,:)); %calculate the mean of the number of neighbors in each class
end

%neighborless targets (targets with no neighbor into the range)
target_idx_zero_neigh_opt = [];
for target = 1 : length(z_cal) %for each target
    if her.n_pairs_by_class_obs(end,target) == sum(her.n_pairs_by_class_obs(2:end,target))
        target_idx_zero_neigh_opt = [target_idx_zero_neigh_opt target];
    end
end

% step 5: z PMF
%create a z PMF for each neighbor based on the neighbor class
pmf_z_target_given_neigh_opt = cell(length(z_cal),length(z_cal)); %z PMF for each neighbor (row) (columns represent analyzed target)
diff_n_bins_shift = her.n_bins_z*her.n_bins_z_per_bin_shift - her.n_bins_shift;

for target = 1 : length(z_cal) %for each target
    for i = 1 : length(z_cal) %for each neighbor
        class_ = classes_obs_target_opt(i,target); %take the distance class
        if class_ ~= 0 % in case there is an associated class 
            bins_shift = round((her.edges_diff_z_shift(1) + z_cal(i) - her.edges_z(1)) / her.binwidth_shift);
            idx_bins_shift = reshape( ...
                  [ones(1,bins_shift) 2:her.n_bins_shift+1 ones(1,int16(diff_n_bins_shift) - bins_shift)], ... %fill edges with 1 (zero bin)
                  int16(her.n_bins_z_per_bin_shift), []);
            pmf_ = her.pmf_diff_z_by_class_obs_range_shift(class_, :);
            pmf_z_target_given_neigh_opt{i,target} = sum(pmf_(idx_bins_shift), 1);
            pmf_z_target_given_neigh_opt{i,target} = pmf_z_target_given_neigh_opt{i,target} / sum(pmf_z_target_given_neigh_opt{i,target}); %normalize
        end
    end
end

% step 6: Convex optimization OR
weights_or = nan(her.n_classes_range,1);
for class_ = 1 : her.n_classes_range
    weights_or(class_,1) = (1 / her.bin_centers_distance_classes_range(class_)) / (sum(1./her.bin_centers_distance_classes_range));  % define weights for each distance class (first guess of the optimization)
end
w_OR0 = weights_or(2:end) / weights_or(1); %normalize to 1st element
fobj = @(w_OR) f_DKL_w_OR([1; w_OR(:)], z_cal', classes_obs_target_opt, pmf_z_target_given_neigh_opt, her.edges_z); %to avoid w multiples of 10, set the 1st class weight to 1 
A = -eye(her.n_classes_range-1) + diag(ones(her.n_classes_range-2,1), 1);  % A*x <= b  ==>  w_OR(i+1) <= w_OR(i)
A = A(1:end-1,:);
b = zeros(her.n_classes_range-2, 1);
lb = 1e-6*ones(her.n_classes_range-1,1); %lower bounds
ub = ones(her.n_classes_range-1,1); %upper bounds
w_OR = fmincon(fobj, w_OR0, A, b, [], [], lb, ub, [], optimoptions('fmincon', 'Display','iter-detailed', 'UseParallel',true));
her.best_w_OR = [1; w_OR(:)] ./ sum([1; w_OR(:)]);

[her.DKL_w_OR, weight_obs_or, pmf_OR] = f_DKL_w_OR(her.best_w_OR, z_cal', classes_obs_target_opt, pmf_z_target_given_neigh_opt, her.edges_z); % step 5: Normalized weights & step 7: OR PMF - Cross-validation

% Create a continuous OR (based on slope) and the weights OR by class
% The slope of the w_OR is used to avoid discrete OR weights
% slope of the last bin (full dataset) equals to the penultimate
for i = 1:length(her.best_w_OR)-1
    her.w_OR_slope(i) = (her.best_w_OR(i+1) - her.best_w_OR(i)) / (her.edges_distance_classes_range(i+1) - her.edges_distance_classes_range(i));    
end
her.w_OR_slope(her.n_classes_range) = 0;

% step 7: Convex optimization AND
%w_OR does not matter here
weights_and = nan(her.n_classes_range,1);
for class_ = 1 : her.n_classes_range
    weights_and(class_,1) = (1 / her.bin_centers_distance_classes_range(class_)) / (sum(1./her.bin_centers_distance_classes_range));  % define weights for each distance class (first guess of the optimization)
end

w_AND0 = weights_and/max(weights_and); %normalize to last element
fobj = @(w_AND) f_DKL_w_AND(w_AND, z_cal', classes_obs_target_opt, pmf_z_target_given_neigh_opt, her.edges_z); %to avoid w multiples of 10, set the 1st class weight to 1 
A = -eye(her.n_classes_range) + diag(ones(her.n_classes_range-1,1), 1);  % A*x <= b  ==>  w(i+1) <= w(i)
A = A(1:end-1,:);
b = zeros(her.n_classes_range-1, 1);
lb = zeros(her.n_classes_range,1); %lower bounds
ub = ones(her.n_classes_range,1) + eps; %upper bounds
w_AND = fmincon(fobj, w_AND0, A, b, [], [], lb, ub, [], optimoptions('fmincon', 'Display','iter-detailed', 'UseParallel',true));
her.best_w_AND = w_AND;

[her.DKL_w_AND, weight_obs_and, pmf_AND] = f_DKL_w_AND(her.best_w_AND, z_cal', classes_obs_target_opt, pmf_z_target_given_neigh_opt, her.edges_z); % step 5: Normalized weights & step 7: OR PMF - Cross-validation

%Create a continuous AND (based on slope) and weights AND by class
% The slope of the w_AND is used to avoid discrete AND weights
% slope of the last bin (full dataset) equals to the penultimate
for i = 1:length(her.best_w_AND)-1
    her.w_AND_slope(i) = (her.best_w_AND(i+1) - her.best_w_AND(i)) / (her.edges_distance_classes_range(i+1) - her.edges_distance_classes_range(i));    
end
her.w_AND_slope(her.n_classes_range) = 0;%her.w_OR_slope(n_classes-1);

% step 8: Grid search Alpha and Beta 
% grid search for the log-linear combination of AND & OR PMFs
alpha = [1:-.05:0]; %create a range of weights for the AND_OR PMF combination (1 means purely AND combination (intersection of PMFs), 0 purely OR)
beta = alpha;
pmf_alpha_beta_nn = cell(length(alpha),length(z_cal)); %cell for the predicted z PMFs (row=set of weights, column=target)

% finding the closest ones
her.n_neighbor_aggreg = her.n_neighbor - 1;

her.idx_opt_nn = NaN(her.n_neighbor_aggreg,length(z_cal));
mat_euc_distance_xy_nn = NaN(her.n_neighbor_aggreg,length(z_cal));
classes_obs_nn = NaN(her.n_neighbor_aggreg,length(z_cal));
z_target_opt_nn = NaN(her.n_neighbor_aggreg,length(z_cal));
pmf_diff_z_plus_z_nn = cell(her.n_neighbor_aggreg,length(z_cal)); 

for target = 1 : length(z_cal) %for each target identidy the closest neighbors
    [~,idx_] = sort(her.mat_euc_distance_xy(:,target));
    her.idx_opt_nn(:, target) = idx_(2:her.n_neighbor_aggreg+1); %avoinding the class 0 (point in relation to itself)
    mat_euc_distance_xy_nn(:,target) = her.mat_euc_distance_xy(her.idx_opt_nn(:, target),target);
    classes_obs_nn(:,target) = classes_obs_target_opt(her.idx_opt_nn(:, target),target);
    z_target_opt_nn(:,target) = z_cal(her.idx_opt_nn(:, target));
    pmf_diff_z_plus_z_nn(:,target) = pmf_z_target_given_neigh_opt(her.idx_opt_nn(:, target),target);     
end

weights_cont_or_obs_nn = zeros(her.n_neighbor_aggreg,length(z_cal));
weights_cont_and_obs_nn = zeros(her.n_neighbor_aggreg,length(z_cal)); 
normalized_weight_or_cont_obs_nn = zeros(her.n_neighbor_aggreg,length(z_cal));

for target = 1:length(z_cal)
    for obs = 1:her.n_neighbor_aggreg
        weights_cont_or_obs_nn(obs,target) =  her.best_w_OR(classes_obs_nn(obs,target)) + her.w_OR_slope(classes_obs_nn(obs,target)) * ( mat_euc_distance_xy_nn(obs,target) - her.edges_distance_classes_range(classes_obs_nn(obs,target)) ); %continuous OR weights
        weights_cont_and_obs_nn(obs,target) =  her.best_w_AND(classes_obs_nn(obs,target)) + her.w_AND_slope(classes_obs_nn(obs,target)) * ( mat_euc_distance_xy_nn(obs,target) - her.edges_distance_classes_range(classes_obs_nn(obs,target)) ); %continuous OR weights

    end
end

for target = 1:length(z_cal)
    for obs = 1:her.n_neighbor_aggreg
        normalized_weight_or_cont_obs_nn(obs,target) = weights_cont_or_obs_nn(obs,target) ./ sum(weights_cont_or_obs_nn(:,target)); %continuous OR weights
     end
end

% calculate alpha_beta PMF for target 
w_alpha_beta = NaN(length(alpha)*length(beta),2);
i=1; %counter
for alpha_ = alpha  %for each set of alpha weights 
    for beta_ = beta %for each set of alpha weights 
        for target = 1 : length(z_cal) %for each target
            idx = [1:her.n_neighbor_aggreg]; %it ignores when i=j
            pmfs_ = cell2mat(pmf_diff_z_plus_z_nn(idx,target)); %temporary vector of the PMF 
            weights_or_ = normalized_weight_or_cont_obs_nn(idx,target); %take the best set of OR weight
            weights_and_ = weights_cont_and_obs_nn(idx,target);
            pmfs_and_ = f_loglinear_aggregation(pmfs_,weights_and_);
            pmfs_or_ = f_linear_aggregation(pmfs_,weights_or_);
            [ pmf_alpha_beta_nn{i,target} ] = f_loglinear_aggregation([pmfs_and_; pmfs_or_], [alpha_; beta_]); %save the predicted z PMF of the target (0 = purely OR combination)
        end
        w_alpha_beta(i,1:2) = [alpha_ beta_];
        i=i+1;
    end
end

probab_z_obs_ = NaN(1,length(z_cal)); %temporary variable
DKL_w_alpha_beta_ = NaN(size(w_alpha_beta,1),1); %DKL between the true Z and the predicted one. Each row corresponds to a set of weights
probab_AND_true_obs = ones(size(w_alpha_beta,1),length(z_cal));

for w = 1 : size(w_alpha_beta,1) %for each set of AND_OR weights
    prob_ = pmf_alpha_beta_nn(w,:);
    PMF_true = ones(1,size(z_cal,1));
    DKL_w_alpha_beta_(w) = f_performance_prob(z_cal', prob_, PMF_true, her.edges_z);
end

idx_best_alpha_beta = min(find(DKL_w_alpha_beta_==min(DKL_w_alpha_beta_))); %em caso de empate, vamos pegar o que tem mais OR
her.DKL_w_alpha_beta = DKL_w_alpha_beta_(idx_best_alpha_beta);
her.best_alpha = w_alpha_beta(idx_best_alpha_beta,1);
her.best_beta = w_alpha_beta(idx_best_alpha_beta,2);
pmf_z_obs_AND_OR = pmf_alpha_beta_nn(idx_best_alpha_beta,:); %the best z PMFs 

her.bin_centers_edges_z = [];
her.bin_centers_edges_z = her.edges_z(1:length(her.edges_z)-1) + her.binwidth_z/2;

end