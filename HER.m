%% Load dataset 
clear all;
close all;
clc; 
dim_cal = 200
txt = 'LR1'
load(sprintf('datasets/%s_input_randomfield_cal%i', txt, dim_cal));
addpath('functions/');
her.txt = txt; % dataset name

%% Define infogram and Geo3 properties
her.lag_dist = 2.0; % DEFINE the class size (lag width in meters if the input x,y is given in meters)
her.binwidth = 0.20; % DEFINE delta z binwidth
her.binwidth_z = her.binwidth; % DEFINE binwidth for the z PMFs
her.n_neighbor = 12; % DEFINE the maximum number of neighbord to be considered in the interpolation and AND_OR optimization (n_neighbor-1) (min 1, max dim_cal)
%to be used in Geo3 part (no need to run Geo1 and Geo2 parts when you change it)
    her.aggregation_type = 'andor';  %DEFINE the aggregation method
    z_thresh = NaN; %DEFINE z threshold for obtaining the probability map 
                    %(probability of exceeding z_thresh threshold). 
                    %NaN if no one is desired

%% Geo1: Spatial Characterization
% step 1: get diff_z as a function of the distance 
% calculate the euclidean distance and the difference between all point pairs
zm_res_ = min(diff(unique(z_cal))); % smallest Z resolution
her.n_bins_z_per_bin_shift = min(100, her.binwidth_z / zm_res_); % limit to 100 or binwidth_shift = 1 (in the case where all zm are integers)
her.binwidth_shift = her.binwidth_z / her.n_bins_z_per_bin_shift;

her.mat_euc_distance_xy = NaN(length(x_cal),length(x_cal)); %matrix for the distance between 2 pair points
her.mat_diff_z = NaN(length(z_cal),length(z_cal)); %matrix for the diff_z 
for i = 1:length(x_cal) %repeat to each east coordinates (x)
    for j = 1:length(x_cal) %repeat to each north coordinates (y)
        her.mat_euc_distance_xy(i,j) = f_euclidean_dist(x_cal(i), y_cal(i), x_cal(j), y_cal(j)); %calculate the euclidean distance
        her.mat_diff_z(i,j) = f_diff(z_cal(i), z_cal(j)); %calculate the diff_z
    end
end

% step 2: define the distance classes (lag)
her.n_lag = ceil( max(her.mat_euc_distance_xy(:)) / her.lag_dist ); %total number of classes
her.edges_distance_classes = [0:her.lag_dist:her.lag_dist*her.n_lag]; %edges of the distance classes
her.bin_centers_distance_classes = her.edges_distance_classes(1:end-1) + her.lag_dist/2; %bin centers of the distance classes

% check empty classes
count = 1;
her.emptyclass = [];
for i = 1 : her.n_lag %for each class
    idx = her.mat_euc_distance_xy > her.edges_distance_classes(i) & her.mat_euc_distance_xy <= her.edges_distance_classes(i+1); %find the diff_z within the current lag class
    if sum(idx(:)) == 0
        her.emptyclass(count) = i;  
        count = count + 1;
    end
end

if ~isempty(her.emptyclass) %if we have empty bins, inform
    str = ['Atention! There are empty distance classes (', num2str(her.emptyclass), '). Consider increase the lag or bin width to fill the classes into the range.'];
    disp(str)
end

% step 3: compute bin edges
max_diff_z = max(her.mat_diff_z(:));
n_bins = 2*floor(max_diff_z/her.binwidth + 1); % symmetric
mini = -n_bins/2*her.binwidth;
maxi =  n_bins/2*her.binwidth;
her.edges_diff_z = linspace(mini,maxi,n_bins+1); %edges of the diff_z pmf

her.n_bins_shift = 2*floor(max_diff_z/her.binwidth_shift + 1); % symmetric
mini = -her.n_bins_shift/2*her.binwidth_shift;
maxi =  her.n_bins_shift/2*her.binwidth_shift;
her.edges_diff_z_shift = linspace(mini,maxi,her.n_bins_shift+1); %edges of the smaller diff_z pmf

% step 4: Compute delta_z PMF for the full dataset
idx = her.mat_euc_distance_xy > 0; %ignore the index of the observation in relation to itself
her.obs_diff_z = her.mat_diff_z(idx); %take the diff_z observations
[pmf_diff_z_all_obs_withemptybins,~] = histcounts(her.obs_diff_z,her.edges_diff_z,'Normalization', 'probability');
[countpmf_diff_z_all_obs,~] = histcounts(her.obs_diff_z,her.edges_diff_z); %number of pair in each bin
countpmf_diff_z_all_obs_plus1 = countpmf_diff_z_all_obs + 1; %avoiding empty bins (PLUS 1 in each bin)
pmf_diff_z_all_obs = countpmf_diff_z_all_obs_plus1 / sum(countpmf_diff_z_all_obs_plus1);
prob_min = 1/(numel(her.obs_diff_z));
% max(pmf_diff_z_all_obs_withemptybins(:) - pmf_diff_z_all_obs(:)) %checking the maximum difference between the true PMF and PMF_plus 1

[countpmf_diff_z_all_obs_shift,~] = histcounts(her.obs_diff_z,her.edges_diff_z_shift);
countpmf_diff_z_all_obs_shift_plus1 = countpmf_diff_z_all_obs_shift + 1/her.n_bins_z_per_bin_shift; %avoiding empty bins (PLUS 1/n in each bin)
countpmf_diff_z_all_obs_shift_plus1 = [1/her.n_bins_z_per_bin_shift countpmf_diff_z_all_obs_shift_plus1]; %add empty bin in beginning
pmf_diff_z_all_obs_shift = countpmf_diff_z_all_obs_shift_plus1 / sum(countpmf_diff_z_all_obs_shift_plus1);
prob_min_shift = prob_min/her.n_bins_z_per_bin_shift;

% check probabilities sum to 1
if abs(sum(pmf_diff_z_all_obs) - 1) > .00001
    error('Probablities dont sum to 1.')
end
if abs(sum(pmf_diff_z_all_obs_shift) - 1) > .00001
    error('Probablities dont sum to 1.')
end

% step 5: Calculate the delta_z PMF by distance class
%address delta_z observations to its distance classes
obs_diff_z_by_class = cell(1, length(her.edges_distance_classes)-1); %observations of diff_z by lag class (each cell is one lag class) 
pmf_diff_z_by_class_obs_withemptybins = nan(length(her.edges_distance_classes)-1, n_bins);
pmf_diff_z_by_class_obs = nan(length(her.edges_distance_classes)-1, n_bins); %diff_z PMF by lag class (each row is one class, the columns are the probability of the diff_z bins)
pmf_diff_z_by_class_obs_shift_withemptybins = nan(length(her.edges_distance_classes)-1, her.n_bins_shift);
pmf_diff_z_by_class_obs_shift = nan(length(her.edges_distance_classes)-1, her.n_bins_shift+1); %diff_z PMF by lag class (each row is one class, the columns are the probability of the diff_z bins)

for i = 1 : length(her.edges_distance_classes)-1 %for each class
    idx = her.mat_euc_distance_xy > her.edges_distance_classes(i) & her.mat_euc_distance_xy <= her.edges_distance_classes(i+1); %find the diff_z within the current lag class
    obs_diff_z_by_class{i} = her.mat_diff_z(idx); %save the diff_z 
    [pmf_diff_z_by_class_obs_withemptybins(i,:),~] = histcounts(obs_diff_z_by_class{i},her.edges_diff_z,'Normalization', 'probability');
    pmf_diff_z_by_class_obs(i,:) = pmf_diff_z_by_class_obs_withemptybins(i,:) + prob_min;
    pmf_diff_z_by_class_obs(i,:) = pmf_diff_z_by_class_obs(i,:) ./ sum(pmf_diff_z_by_class_obs(i,:)) ; %normalization

    [pmf_diff_z_by_class_obs_shift_withemptybins(i,:),~] = histcounts(obs_diff_z_by_class{i},her.edges_diff_z_shift,'Normalization', 'probability');
    pmf_diff_z_by_class_obs_shift(i,2:end) = pmf_diff_z_by_class_obs_shift_withemptybins(i,:) + prob_min_shift;
    pmf_diff_z_by_class_obs_shift(i,2:end) = pmf_diff_z_by_class_obs_shift(i,2:end) ./ sum(pmf_diff_z_by_class_obs_shift(i,2:end)) ;
end
pmf_diff_z_by_class_obs_shift(:,1) = prob_min_shift; %add empty bin in beginning
% max(max(pmf_diff_z_by_class_obs(1:end-1,:) - pmf_diff_z_by_class_obs_withemptybins(1:end-1,:))) %check maximum diff of the true PMF and PMF_plus 1

% check probabilities sum to 1
for i = 1 : length(her.edges_distance_classes)-1 %for each class
    if abs(sum(pmf_diff_z_by_class_obs(i,:)) - 1) > .00001
        error('Probablities dont sum to 1.')
    end
    if abs(sum(pmf_diff_z_by_class_obs_shift(i,2:end)) - 1) > .00001
        error('Probablities dont sum to 1.')
    end
end

% step 6: Entropy
%calculate the entropy of the delta_z PMFs for the full dataset and for the distance classes 
her.H_diff_z = f_entropy(pmf_diff_z_all_obs); %entropy of the diff_z PMF of the full dataset
her.H_diff_z_by_class = NaN(length(her.edges_distance_classes)-1,1); %vector for the entropy of the diff_z PMF by lag classes 

for i = 1:length(her.edges_distance_classes)-1 %for each lag class
    her.H_diff_z_by_class(i,1) = f_entropy(pmf_diff_z_by_class_obs(i,:)); %calculate the entropy
end

% step 7: Define the info RANGE, and associate delta_z PMF of the full dataset to the distance classes beyond the range
%calculate the PMF of the full dataset where the class entropy > full dataset entropy 
%info RANGE: limit where the entropy of the lag class is greater than the entropy of the full dataset
if her.H_diff_z > her.H_diff_z_by_class(~isnan(her.H_diff_z_by_class))
    her.n_classes_limit = 15;  %DEFINE it looking at the infogram
    her.H_lim = her.H_diff_z_by_class(her.n_classes_limit);
    her.n_classes_range =  find(her.H_diff_z_by_class > her.H_lim,1);
else
    her.H_lim = her.H_diff_z;
    her.n_classes_range = find(her.H_diff_z_by_class > her.H_lim,1); %number of distance classes inside the range + full dataset
end

her.edges_distance_classes_range = [0 ]; %edges of distance classes considering the range (starting in zero meters)
her.bin_centers_distance_classes_range = [];
her.pmf_diff_z_by_class_obs_range = []; %diff_z PMF by lag class considering the range (each row is one class, the columns are the probability of the diff_z bins)
her.H_diff_z_by_class_range = []; %entropy of the delta_z PMF by lag classes inside the range + full dataset
obs_diff_z_by_class_range = []; %observations of diff_z by lag class considering the range (each cell is one lag class) 

for i = 1 : size(pmf_diff_z_by_class_obs,1) %for each class
    if her.H_diff_z_by_class(i,1) > her.H_lim %stop if the entropy of the lag class is greater than the entropy of the full dataset (or a defined one)
        break
    else
        her.pmf_diff_z_by_class_obs_range(i,:) = pmf_diff_z_by_class_obs(i,:); %save the diff_z PMF
        her.pmf_diff_z_by_class_obs_range_shift(i,:) = pmf_diff_z_by_class_obs_shift(i,:); %save the diff_z PMF
        her.edges_distance_classes_range(1,i+1) = her.edges_distance_classes(1,i+1); %save the edge of the distance class   
        her.H_diff_z_by_class_range(1,i) = her.H_diff_z_by_class(i); %save the entropy of the diff_z PMF
        obs_diff_z_by_class_range{i} = obs_diff_z_by_class{i}; %save the diff_z observation
    end
end

%create a last distance class and associate to it 
her.pmf_diff_z_by_class_obs_range(size(her.pmf_diff_z_by_class_obs_range,1)+1,:) = pmf_diff_z_all_obs; %diff_z PMF of the full dataset
her.pmf_diff_z_by_class_obs_range_shift(size(her.pmf_diff_z_by_class_obs_range_shift,1)+1,:) = pmf_diff_z_all_obs_shift; %diff_z PMF of the full dataset
her.edges_distance_classes_range(1,length(her.edges_distance_classes_range)+1) = her.edges_distance_classes(1,end); %the right edge of the distance class 
her.H_diff_z_by_class_range(1,length(her.H_diff_z_by_class_range)+1) = her.H_diff_z; %the entropy of the diff_z PMF of the full dataset
obs_diff_z_by_class_range{length(obs_diff_z_by_class_range)+1} = her.obs_diff_z; %the diff_z of the full dataset

%calculate bin centers of the distance classes
for i = 1:length(her.edges_distance_classes_range)-1
    her.bin_centers_distance_classes_range(i) = 1/2 * (her.edges_distance_classes_range(i) + her.edges_distance_classes_range(i+1));
end

%distribution of the # of pairs by class (used for the histogram construction)
for i = 1:her.n_lag
    her.n_pairs_by_class(i) = numel(obs_diff_z_by_class{i});
end

if isnan(sum(her.pmf_diff_z_by_class_obs_range(:)))
    str = ['Atention! There are no histograms for \Deltaz in the range (classes without histogram: ' num2str(find( sum(her.pmf_diff_z_by_class_obs_range(:,:)') ~= 1 )) '). Interpolation not possivel. Consider increase \Deltaz binwidth.'];
    disp(str)
end  

% plot Infogram cloud, Infogram & Number of pairs, Histogram of delta_z classes, normalized infogram + variogram
f_plot_infogram(her);

%% Geo2: Weight optimization 
% (leave-one-out cross validation on calibration dataset)
% step 1: define edges of Z 
max_z = max(z_cal(:));
min_z = min(z_cal(:));
max_diff_z = max(her.obs_diff_z(:));
min_diff_z = min(her.obs_diff_z(:));

% compute bin edges
n_bins_left = -floor((min_diff_z + min_z)/ her.binwidth_z);
n_bins_right = floor((max_diff_z + max_z)/ her.binwidth_z + 1); %edge <= x < edge
her.n_bins_z = n_bins_left+n_bins_right;
mini = -her.binwidth_z*n_bins_left;
maxi =  her.binwidth_z*n_bins_right;
her.edges_z = linspace(mini,maxi,her.n_bins_z+1); %edges of the z+diff_z pmf
[pmf_z_all_obs,~] = histcounts(z_cal, her.edges_z,'Normalization', 'probability'); %calculate the z PMF of the full dataset 

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
            pmf_z_target_given_neigh_opt{i,target} = sum(pmf_(idx_bins_shift));
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

% plot optimum weight
f_plot_weights(her);

%% Geo3: z PMF prediction
%GRID for ploting results
[x_target_grid,y_target_grid] = meshgrid(1:100, 1:100); %DEFINE a grid for predicting and ploting results
x_target = x_target_grid(:);
y_target = y_target_grid(:);

%Predict GRID
[pmf_pred_nn, target_idx_zero_neigh_pred] = f_her_predict(x_cal, y_cal,z_cal, x_target, y_target, her);

% %Predict validation set (for performance analysis, if your validation set does not match with the grid)
% % x and y here could be a column vector
% [pmf_pred_nn, target_idx_zero_neigh_pred] = f_her_predict(x_cal, y_cal,z_cal, x, y, her);

%% Extract PMF statistics
% extract mean, median, mode, probability and plot

% GRID: plot maps and predicted PMFs
[z_target_entropy_pred_GRID_plot, z_target_mean_pred_GRID_plot, z_target_median_pred_GRID_plot, z_target_mode_pred_GRID_plot, z_target_probability_pred_GRID_plot] = f_extract_pmf_statistics(x_target_grid, y_target_grid, pmf_pred_nn, her.bin_centers_edges_z, z_thresh);
x_pmf = [10, 25, 47, 49]; y_pmf = [42, 63, 16, 73]; %coordinates to plot predicted PMF
f_plot_prediction(z_target_mean_pred_GRID_plot, z_target_entropy_pred_GRID_plot, pmf_pred_nn, x_target, y_target, x_target_grid, y_target_grid, x, y, z, idx_cal, her, x_pmf, y_pmf,'Y');
% f_plot_probabilitymap(z_target_probability_pred_GRID_plot, z_thresh, txt, x_target_grid, y_target_grid, x, y, z, idx_cal);

%% Calculate performance metrics
% Root mean square error (RMSE)
% Mean Error (ME)  
% Mean Absolute Error (MAE) L1 norm, robust parameter estimator
% Nash-Sutcliffe model efficiency (r2, coefficient of determination)
z_target_pred_ = z_target_mean_pred_GRID_plot;
pmf_ = pmf_pred_nn;

% performance validation set
[perf.error_sign_val, perf.RMSE_val, perf.ME_val, perf.MAE_val, perf.NSE_val] = f_performance_det(z_target_pred_(idx_val),z(idx_val));
% scoring rule - DKL
PMF_true_val = ones(1,length(idx_val));
PMF_simulated_val = pmf_(idx_val);
perf.DKL_score_mean_val = f_performance_prob(z_val', PMF_simulated_val, PMF_true_val, her.edges_z);
perf.correl_val = corr(z_val, (z_target_pred_(idx_val) - z_val));

%%% performance test set
[perf.error_sign_test, perf.RMSE_test, perf.ME_test, perf.MAE_test, perf.NSE_test] = f_performance_det(z_target_pred_(idx_test),z(idx_test));
% scoring rule - DKL
PMF_true_test = ones(1,length(idx_test));
PMF_simulated_test = pmf_(idx_test);
perf.DKL_score_mean_test = f_performance_prob(z_test', PMF_simulated_test, PMF_true_test, her.edges_z);
perf.correl_test = corr(z_test, (z_target_pred_(idx_test) - z_test));

%% clear 
%Geo1
clear count str i j idx im jm k center_left center_right max_diff_z ...
    min_diff_z maxi mini ans left_color right_color samp ...
    pmf_diff_z_all_obs_withemptybins pmf_diff_z_by_class_obs_shift_withemptybins ...
    pmf_diff_z_by_class_obs_withemptybins pmf_diff_z_by_class_obs zm_res_ ...
    pmf_diff_z_all_obs_shift pmf_diff_z_all_obs pmf_diff_z_by_class_obs_shift ...
    obs_diff_z_by_class_range obs_diff_z_by_class n_bins countpmf_diff_z_all_obs_shift_plus1...
    countpmf_diff_z_all_obs countpmf_diff_z_all_obs_plus1 countpmf_diff_z_all_obs_shift   
%Geo2
clear pmf_alpha_beta_nn alpha_ beta_ alpha beta prob_min prob_min_shift ...
    left_color right_color i bins_shift class_ diff_n_bins_shift classes_obs_target_opt...
    classes_obs_nn weights_or_ weights_and_ w DKL_w_alpha_beta_ w_AND0 w_OR0 w_OR w_AND ...
    w_alpha_beta weights_and weights_or weights_cont_and_obs_nn weights_cont_or_obs_nn...
    normalized_weight_or_cont_obs_nn prob_ pmfs_or_ pmfs_and_ pmfs_ probab_z_obs_ pmf_...
    str z_target_opt_nn n_bins_left n_bins_right mini maxi min_diff_z max_z min_z ...
    idx_best_alpha_beta idx_bins_shift mat_euc_distance_obs_target_opt ...
    mat_euc_distance_xy_nn max_diff_z pmf_AND pmf_OR pmf_diff_z_plus_z_nn PMF_true...
    pmf_z_all_obs pmf_z_obs_AND_OR pmf_z_target_given_neigh_opt pmf_z_target_given_neigh_opt...
    target ub lb fobj A b obs ax fig idx idx_ ncols nrows probab_AND_true_obs...
    target_idx_zero_neigh_opt weight_obs_or weight_obs_and z_target_pred_
%Geo3
clear PMF_simulated_val PMF_true_val
