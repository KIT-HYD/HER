clear all;
close all;
clc; 

%% Load dataset
dim_cal = 200
txt = 'LR1'
load(sprintf('datasets/%s_input_randomfield_cal%i', txt, dim_cal));

lag_dist = 2.0; % DEFINE the class size (lag width in meters if the input x,y is given in meters)
binwidth = 0.20; % DEFINE delta z binwidth
binwidth_z = binwidth; % DEFINE binwidth for the z PMFs
n_neighbor = 12; % define the maximum number of neighbord to be considered in the interpolation and AND_OR optimization (n_neighbor-1) (min 1, max dim_cal)

zm_res_ = min(diff(unique(z_cal))); % smallest Z resolution
n_bins_z_per_bin_shift = min(100, binwidth_z / zm_res_); % limit to 100 or binwidth_shift = 1 (in the case where all zm are integers)
binwidth_shift = binwidth_z / n_bins_z_per_bin_shift;

%% Geo1 Spatial Characterization
% step 1: get diff_z as a function of the distance 
% calculate the euclidean distance and the difference between all point pairs
mat_euc_distance_xy = NaN(length(x_cal),length(x_cal)); %matrix for the distance between 2 pair points
mat_diff_z = NaN(length(z_cal),length(z_cal)); %matrix for the diff_z 

for i = 1:length(x_cal) %repeat to each east coordinates (x)
    for j = 1:length(x_cal) %repeat to each north coordinates (y)
        mat_euc_distance_xy(i,j) = f_euclidean_dist(x_cal(i), y_cal(i), x_cal(j), y_cal(j)); %calculate the euclidean distance
        mat_diff_z(i,j) = f_diff(z_cal(i), z_cal(j)); %calculate the diff_z
    end
end

% step 2: define the distance classes (lag)
n_lag = ceil( max(mat_euc_distance_xy(:)) / lag_dist ); %total number of classes
edges_distance_classes = [0:lag_dist:lag_dist*n_lag]; %edges of the distance classes
bin_centers_distance_classes = edges_distance_classes(1:end-1) + lag_dist/2; %bin centers of the distance classes

% check empty classes
count = 1;
emptyclass = [];
for i = 1 : n_lag %for each class
    idx = mat_euc_distance_xy > edges_distance_classes(i) & mat_euc_distance_xy <= edges_distance_classes(i+1); %find the diff_z within the current lag class
    if sum(idx(:)) == 0
        emptyclass(count) = i;  
        count = count + 1;
    end
end

if ~isempty(emptyclass) %if we have empty bins, inform
    str = ['Atention! There are empty distance classes (', num2str(emptyclass), '). Consider increase the lag or bin width to fill the classes into the range.'];
    disp(str)
end

% EXTRA step: create a experimental semivariogram
lag_variance = NaN(length(edges_distance_classes)-1,1); %vector for the semivariance of the observations
n_obs_by_lag = zeros(length(edges_distance_classes)-1,1); %vector for counting the number of pair points in each lag class

mat_square_diff_z = mat_diff_z.^2; %matrix of the squared semivariance

for i = 1 : (length(edges_distance_classes)-1) %for each distance class
    idx = mat_euc_distance_xy > edges_distance_classes(i) & mat_euc_distance_xy <= edges_distance_classes(i+1); %find the observations within the current lag class
    n_obs_by_lag(i)= sum(idx(:)); %count the number of values
    lag_variance(i) = 1/2 * mean(mat_square_diff_z(idx)); %calculate the semivariance
end

% step 3: compute bin edges
max_diff_z = max(mat_diff_z(:));
n_bins = 2*floor(max_diff_z/binwidth + 1); % symmetric
mini = -n_bins/2*binwidth;
maxi =  n_bins/2*binwidth;
edges_diff_z = linspace(mini,maxi,n_bins+1); %edges of the diff_z pmf

n_bins_shift = 2*floor(max_diff_z/binwidth_shift + 1); % symmetric
mini = -n_bins_shift/2*binwidth_shift;
maxi =  n_bins_shift/2*binwidth_shift;
edges_diff_z_shift = linspace(mini,maxi,n_bins_shift+1); %edges of the smaller diff_z pmf

% step 4: Compute delta_z PMF for the full dataset
idx = mat_euc_distance_xy > 0; %ignore the index of the observation in relation to itself
obs_diff_z = mat_diff_z(idx); %take the diff_z observations
[pmf_diff_z_all_obs_withemptybins,~] = histcounts(obs_diff_z,edges_diff_z,'Normalization', 'probability');
[countpmf_diff_z_all_obs,~] = histcounts(obs_diff_z,edges_diff_z); %number of pair in each bin
countpmf_diff_z_all_obs_plus1 = countpmf_diff_z_all_obs + 1; %avoiding empty bins (PLUS 1 in each bin)
pmf_diff_z_all_obs = countpmf_diff_z_all_obs_plus1 / sum(countpmf_diff_z_all_obs_plus1);
prob_min = 1/(numel(obs_diff_z));
max(pmf_diff_z_all_obs_withemptybins(:) - pmf_diff_z_all_obs(:)) %checking the maximum difference between the true PMF and PMF_plus 1

[countpmf_diff_z_all_obs_shift,~] = histcounts(obs_diff_z,edges_diff_z_shift);
countpmf_diff_z_all_obs_shift_plus1 = countpmf_diff_z_all_obs_shift + 1/n_bins_z_per_bin_shift; %avoiding empty bins (PLUS 1/n in each bin)
countpmf_diff_z_all_obs_shift_plus1 = [1/n_bins_z_per_bin_shift countpmf_diff_z_all_obs_shift_plus1]; %add empty bin in beginning
pmf_diff_z_all_obs_shift = countpmf_diff_z_all_obs_shift_plus1 / sum(countpmf_diff_z_all_obs_shift_plus1);
prob_min_shift = prob_min/n_bins_z_per_bin_shift;

% check probabilities sum to 1
if abs(sum(pmf_diff_z_all_obs) - 1) > .00001
    error('Probablities dont sum to 1.')
end
if abs(sum(pmf_diff_z_all_obs_shift) - 1) > .00001
    error('Probablities dont sum to 1.')
end

% step 5: Calculate the delta_z PMF by distance class
%address delta_z observations to its distance classes
obs_diff_z_by_class = cell(1, length(edges_distance_classes)-1); %observations of diff_z by lag class (each cell is one lag class) 
pmf_diff_z_by_class_obs_withemptybins = nan(length(edges_distance_classes)-1, n_bins);
pmf_diff_z_by_class_obs = nan(length(edges_distance_classes)-1, n_bins); %diff_z PMF by lag class (each row is one class, the columns are the probability of the diff_z bins)
pmf_diff_z_by_class_obs_shift_withemptybins = nan(length(edges_distance_classes)-1, n_bins_shift);
pmf_diff_z_by_class_obs_shift = nan(length(edges_distance_classes)-1, n_bins_shift+1); %diff_z PMF by lag class (each row is one class, the columns are the probability of the diff_z bins)

for i = 1 : length(edges_distance_classes)-1 %for each class
    idx = mat_euc_distance_xy > edges_distance_classes(i) & mat_euc_distance_xy <= edges_distance_classes(i+1); %find the diff_z within the current lag class
    obs_diff_z_by_class{i} = mat_diff_z(idx); %save the diff_z 
    [pmf_diff_z_by_class_obs_withemptybins(i,:),~] = histcounts(obs_diff_z_by_class{i},edges_diff_z,'Normalization', 'probability');
    pmf_diff_z_by_class_obs(i,:) = pmf_diff_z_by_class_obs_withemptybins(i,:) + prob_min;
    pmf_diff_z_by_class_obs(i,:) = pmf_diff_z_by_class_obs(i,:) ./ sum(pmf_diff_z_by_class_obs(i,:)) ; %normalization

    [pmf_diff_z_by_class_obs_shift_withemptybins(i,:),~] = histcounts(obs_diff_z_by_class{i},edges_diff_z_shift,'Normalization', 'probability');
    pmf_diff_z_by_class_obs_shift(i,2:end) = pmf_diff_z_by_class_obs_shift_withemptybins(i,:) + prob_min_shift;
    pmf_diff_z_by_class_obs_shift(i,2:end) = pmf_diff_z_by_class_obs_shift(i,2:end) ./ sum(pmf_diff_z_by_class_obs_shift(i,2:end)) ;
end
pmf_diff_z_by_class_obs_shift(:,1) = prob_min_shift; %add empty bin in beginning

max(max(pmf_diff_z_by_class_obs(1:end-1,:) - pmf_diff_z_by_class_obs_withemptybins(1:end-1,:))) %check maximum diff of the true PMF and PMF_plus 1
% check probabilities sum to 1
for i = 1 : length(edges_distance_classes)-1 %for each class
    if abs(sum(pmf_diff_z_by_class_obs(i,:)) - 1) > .00001
        error('Probablities dont sum to 1.')
    end
    if abs(sum(pmf_diff_z_by_class_obs_shift(i,2:end)) - 1) > .00001
        error('Probablities dont sum to 1.')
    end
end

% step 6: Entropy
%calculate the entropy of the delta_z PMFs for the full dataset and for the distance classes 
H_diff_z = f_entropy(pmf_diff_z_all_obs); %entropy of the diff_z PMF of the full dataset
H_diff_z_by_class = NaN(length(edges_distance_classes)-1,1); %vector for the entropy of the diff_z PMF by lag classes 

for i = 1:length(edges_distance_classes)-1 %for each lag class
    H_diff_z_by_class(i,1) = f_entropy(pmf_diff_z_by_class_obs(i,:)); %calculate the entropy
end

% plot Infogram + Variogram
fig = figure();
left_color = [0 0 0.4];
right_color = [.7 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

plot(bin_centers_distance_classes, H_diff_z_by_class./max(H_diff_z_by_class),'Marker', '.', 'MarkerSize', 20);

yyaxis left
ylabel('Normalized entropy [bit]');
xlim([0 bin_centers_distance_classes(end-1)]);
ylim([0 1]);

yyaxis right
plot(bin_centers_distance_classes, lag_variance./(max(lag_variance)),'Marker', '.', 'MarkerSize', 20);
xlabel('lag [m]');
ylabel('Normalized dissimilarity (1/2*mean(delta z²))');
hold on 
for i = 2 : length(edges_distance_classes)
    line([edges_distance_classes(i),edges_distance_classes(i)], get(gca, 'ylim'),'Color','black','LineStyle','--');
end
ylim([0 1]);
title('Normalized Infogram and Experimental variogram');

% step 7: Define the info RANGE, and associate delta_z PMF of the full dataset to the distance classes beyond the range
%calculate the PMF of the full dataset where the class entropy > full dataset entropy 
%info RANGE: limit where the entropy of the lag class is greater than the 
%entropy of the full dataset
if H_diff_z > H_diff_z_by_class(~isnan(H_diff_z_by_class))
    n_classes_limit = 15;  %DEFINE it looking at the infogram
    H_lim = H_diff_z_by_class(n_classes_limit);
    n_classes_range =  find(H_diff_z_by_class > H_lim,1);
else
    H_lim = H_diff_z;
    n_classes_range = find(H_diff_z_by_class > H_lim,1); %number of distance classes inside the range + full dataset
end

edges_distance_classes_range = [0 ]; %edges of distance classes considering the range (starting in zero meters)
bin_centers_distance_classes_range = [];
pmf_diff_z_by_class_obs_range = []; %diff_z PMF by lag class considering the range (each row is one class, the columns are the probability of the diff_z bins)
H_diff_z_by_class_range = []; %entropy of the delta_z PMF by lag classes inside the range + full dataset
obs_diff_z_by_class_range = []; %observations of diff_z by lag class considering the range (each cell is one lag class) 

for i = 1 : size(pmf_diff_z_by_class_obs,1) %for each class
    if H_diff_z_by_class(i,1) > H_lim %stop if the entropy of the lag class is greater than the entropy of the full dataset (or a defined one)
        break
    else
        pmf_diff_z_by_class_obs_range(i,:) = pmf_diff_z_by_class_obs(i,:); %save the diff_z PMF
        pmf_diff_z_by_class_obs_range_shift(i,:) = pmf_diff_z_by_class_obs_shift(i,:); %save the diff_z PMF
        edges_distance_classes_range(1,i+1) = edges_distance_classes(1,i+1); %save the edge of the distance class   
        H_diff_z_by_class_range(1,i) = H_diff_z_by_class(i); %save the entropy of the diff_z PMF
        obs_diff_z_by_class_range{i} = obs_diff_z_by_class{i}; %save the diff_z observation
    end
end

%create a last distance class and associate to it 
pmf_diff_z_by_class_obs_range(size(pmf_diff_z_by_class_obs_range,1)+1,:) = pmf_diff_z_all_obs; %diff_z PMF of the full dataset
pmf_diff_z_by_class_obs_range_shift(size(pmf_diff_z_by_class_obs_range_shift,1)+1,:) = pmf_diff_z_all_obs_shift; %diff_z PMF of the full dataset
edges_distance_classes_range(1,length(edges_distance_classes_range)+1) = edges_distance_classes(1,end); %the right edge of the distance class 
H_diff_z_by_class_range(1,length(H_diff_z_by_class_range)+1) = H_diff_z; %the entropy of the diff_z PMF of the full dataset
obs_diff_z_by_class_range{length(obs_diff_z_by_class_range)+1} = obs_diff_z; %the diff_z of the full dataset

%calculate bin centers of the distance classes
for i = 1:length(edges_distance_classes_range)-1
    bin_centers_distance_classes_range(i) = 1/2 * (edges_distance_classes_range(i) + edges_distance_classes_range(i+1));
end


%plot the distribution of the # of pairs by class (used for the histogram construction)
for i = 1:n_lag
    n_pairs_by_class(i) = numel(obs_diff_z_by_class{i});
end

if isnan(sum(pmf_diff_z_by_class_obs_range(:)))
    str = ['Atention! There are no histograms for delta z in the range (classes without histogram: ' num2str(find( sum(pmf_diff_z_by_class_obs_range(:,:)') ~= 1 )) '). Interpolation not possivel. Consider increase delta z binwidth.'];
    disp(str)
end

% clear
clear count str i j idx im jm k center_left center_right max_diff_z ...
    min_diff_z maxi mini ans left_color right_color samp ...
    pmf_diff_z_all_obs_withemptybins pmf_diff_z_by_class_obs_shift_withemptybins ...
    pmf_diff_z_by_class_obs_withemptybins pmf_diff_z_by_class_obs

%% Geo2 Weight optimization (cross validation on the training dataset)
% step 1: define edges of Z 
max_z = max(z_cal(:));
min_z = min(z_cal(:));
max_diff_z = max(obs_diff_z(:));
min_diff_z = min(obs_diff_z(:));

% compute bin edges
n_bins_left = -floor((min_diff_z + min_z)/ binwidth_z);
n_bins_right = floor((max_diff_z + max_z)/ binwidth_z + 1); %edge <= x < edge
n_bins_z = n_bins_left+n_bins_right;
mini = -binwidth_z*n_bins_left;
maxi =  binwidth_z*n_bins_right;
edges_z = linspace(mini,maxi,n_bins_z+1); %edges of the z+diff_z pmf
[pmf_z_all_obs,~] = histcounts(z_cal, edges_z,'Normalization', 'probability'); %calculate the z PMF of the full dataset 

% step 2: calculate the distance of the target to its neighbors 
% calculate the euclidean distance between target and its neighbors
% target: observation to be analized (predicted and validated)
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
        for class_ = 1 : n_classes_range %for each class
            if mat_euc_distance_obs_target_opt(i,target) > edges_distance_classes_range(class_) & mat_euc_distance_obs_target_opt(i,target) <= edges_distance_classes_range(class_+1) % in case it is within the current lag class
                classes_obs_target_opt(i,target) = class_; %save the identified class
            end          
        end
    end    
end

% step 4: neighbor cardinality calculation
%plot the mean of the # of measurement neighbors by class
nmean_obs_by_class = NaN(n_classes_range,1);
[n_obs_by_class, edges_n_obs] = hist(classes_obs_target_opt(:,:), -0.5:1:n_classes_range); %getting the # of neighbors by class (row) for each target (column)
for class_ = 1 : n_classes_range+1 %for each class (into the range + out of the range)
    nmean_obs_by_class(class_) = mean(n_obs_by_class(class_,:)); %calculate the mean of the number of neighbors in each class
end

%plot the distribution of the # of neighbors (for all targets)
edges_n_obs_all = -0.5:1:(max(n_obs_by_class(2:n_classes_range,:)+1,[],'all'));
edges_n_obs_all_centralbin = edges_n_obs_all(1:end-1) + 0.5;
[pmf_n_obs_all] = histcounts(n_obs_by_class(2:n_classes_range,:), edges_n_obs_all,'Normalization', 'probability'); 

%plot the distribution of the # of neighbors by class
edges_n_obs = -0.5:1:(max(n_obs_by_class(2:3,:),[],'all'));
edges_n_obs_centralbin = edges_n_obs(1:end-1) + 0.5;
[pmf_n_obs_by_class1] = histcounts(n_obs_by_class(2,:), edges_n_obs,'Normalization', 'probability'); %1st class is line 2 (line 1 is class 0 = the target in relation to itself)
[pmf_n_obs_by_class2] = histcounts(n_obs_by_class(3,:), edges_n_obs,'Normalization', 'probability');
% [pmf_n_obs_by_class3] = histcounts(n_obs_by_class(4,:), edges_n_obs,'Normalization', 'probability');
data = [pmf_n_obs_by_class1; pmf_n_obs_by_class2];%; pmf_n_obs_by_class3];

%plot the neighborless targets (targets with no neighbor into the range)
target_idx_zero_neigh_opt = [];
for target = 1 : length(z_cal) %for each target
    if n_obs_by_class(end,target) == sum(n_obs_by_class(2:end,target))
        target_idx_zero_neigh_opt = [target_idx_zero_neigh_opt target];
    end
end

% step 5: z PMF
%create a z PMF for each neighbor based on the neighbor class
pmf_z_target_given_neigh_opt = cell(length(z_cal),length(z_cal)); %z PMF for each neighbor (row) (columns represent analyzed target)
diff_n_bins_shift = n_bins_z*n_bins_z_per_bin_shift - n_bins_shift;

for target = 1 : length(z_cal) %for each target
    for i = 1 : length(z_cal) %for each neighbor
        class_ = classes_obs_target_opt(i,target); %take the distance class
        if class_ ~= 0 % in case there is an associated class 
            bins_shift = round((edges_diff_z_shift(1) + z_cal(i) - edges_z(1)) / binwidth_shift);
            idx_bins_shift = reshape( ...
                  [ones(1,bins_shift) 2:n_bins_shift+1 ones(1,int16(diff_n_bins_shift) - bins_shift)], ... %fill edges with 1 (zero bin)
                  int16(n_bins_z_per_bin_shift), []);
            pmf_ = pmf_diff_z_by_class_obs_range_shift(class_, :);
            pmf_z_target_given_neigh_opt{i,target} = sum(pmf_(idx_bins_shift));
            pmf_z_target_given_neigh_opt{i,target} = pmf_z_target_given_neigh_opt{i,target} / sum(pmf_z_target_given_neigh_opt{i,target}); %normalize
        end
    end
end

% step 6: Convex optimization OR
weights_or = nan(n_classes_range,1);
for class_ = 1 : n_classes_range
    weights_or(class_,1) = (1 / bin_centers_distance_classes_range(class_)) / (sum(1./bin_centers_distance_classes_range));  % define weights for each distance class (first guess of the optimization)
end
w_OR0 = weights_or(2:end) / weights_or(1); %normalize to 1st element
fobj = @(w_OR) f_DKL_w_OR([1; w_OR(:)], z_cal', classes_obs_target_opt, pmf_z_target_given_neigh_opt, edges_z); %to avoid w multiples of 10, set the 1st class weight to 1 
A = -eye(n_classes_range-1) + diag(ones(n_classes_range-2,1), 1);  % A*x <= b  ==>  w_OR(i+1) <= w_OR(i)
A = A(1:end-1,:);
b = zeros(n_classes_range-2, 1);
lb = 1e-6*ones(n_classes_range-1,1); %lower bounds
ub = ones(n_classes_range-1,1); %upper bounds
w_OR = fmincon(fobj, w_OR0, A, b, [], [], lb, ub, [], optimoptions('fmincon', 'Display','iter-detailed', 'UseParallel',true));
best_w_OR = [1; w_OR(:)] ./ sum([1; w_OR(:)]);

[DKL_w_OR, normalized_weight_obs, pmf_OR] = f_DKL_w_OR(best_w_OR, z_cal', classes_obs_target_opt, pmf_z_target_given_neigh_opt, edges_z); % step 5: Normalized weights & step 7: OR PMF - Cross-validation

% Create a continuous OR (based on slope) and Plot of the weights OR by class
%The slope of the W_or is used to avoid discrete OR weights
%slope of the last bin (full dataset) equals to the penultimate
for i = 1:length(best_w_OR)-1
    w_OR_slope(i) = (best_w_OR(i+1) - best_w_OR(i)) / (edges_distance_classes_range(i+1) - edges_distance_classes_range(i));    
end
w_OR_slope(n_classes_range) = 0;

% step 7: Convex optimization AND
%w_OR does not matter here, it is just an input for the f_mixpdf
weights_and = nan(n_classes_range,1);
for class_ = 1 : n_classes_range
    weights_and(class_,1) = (1 / bin_centers_distance_classes_range(class_)) / (sum(1./bin_centers_distance_classes_range));  % define weights for each distance class (first guess of the optimization)
end

w_AND0 = weights_and/max(weights_and); %normalize to last element
fobj = @(w_AND) f_DKL_w_AND(w_AND, z_cal', classes_obs_target_opt, pmf_z_target_given_neigh_opt, edges_z); %to avoid w multiples of 10, set the 1st class weight to 1 
A = -eye(n_classes_range) + diag(ones(n_classes_range-1,1), 1);  % A*x <= b  ==>  w(i+1) <= w(i)
A = A(1:end-1,:);
b = zeros(n_classes_range-1, 1);
lb = zeros(n_classes_range,1); %lower bounds
ub = ones(n_classes_range,1) + eps; %upper bounds
w_AND = fmincon(fobj, w_AND0, A, b, [], [], lb, ub, [], optimoptions('fmincon', 'Display','iter-detailed', 'UseParallel',true));
best_w_AND = w_AND;%[w_AND(:)] ./ sum([w_AND(:)]);

[DKL_w_AND, pmf_AND] = f_DKL_w_AND(best_w_AND, z_cal', classes_obs_target_opt, pmf_z_target_given_neigh_opt, edges_z); % step 5: Normalized weights & step 7: OR PMF - Cross-validation

%Create a continuous OR (based on slope) and Plot of the weights OR by class
%The slope pf the W or is used to avoid discrete OR weights
%slope of the last bin (full dataset) equals to the penultimate

for i = 1:length(best_w_AND)-1
    w_AND_slope(i) = (best_w_AND(i+1) - best_w_AND(i)) / (edges_distance_classes_range(i+1) - edges_distance_classes_range(i));    
end

w_AND_slope(n_classes_range) = 0;%w_OR_slope(n_classes-1); %adjusted!

%  step 8: Alpha and Beta PMF 
%grid search for the linear combination of AND & OR combination of PMFs
alpha = [1:-.05:0]; %create a range of weights for the AND_OR PMF combination (1 means purely AND combination (intersection of PMFs), 0 purely OR)
beta = alpha;

pmf_alpha_beta_nn = cell(length(alpha),length(z_cal)); %cell for the predicted z PMFs (row=set of weights, column=target)

% finding the closest ones TODO check it
n_neighbor_andor = n_neighbor - 1;

idx_opt_nn = NaN(n_neighbor_andor,length(z_cal));
mat_euc_distance_xy_nn = NaN(n_neighbor_andor,length(z_cal));
classes_obs_nn = NaN(n_neighbor_andor,length(z_cal));
z_target_opt_nn = NaN(n_neighbor_andor,length(z_cal)); %mesmo do geo 3?
pmf_diff_z_plus_z_nn = cell(n_neighbor_andor,length(z_cal)); 

for target = 1 : length(z_cal) %for each target identidy the closest neighbors
    [~,idx_] = sort(mat_euc_distance_xy(:,target));
    idx_opt_nn(:, target) = idx_(2:n_neighbor_andor+1); %avoinding the class 0 (point in relation to itself)
    mat_euc_distance_xy_nn(:,target) = mat_euc_distance_xy(idx_opt_nn(:, target),target);
    classes_obs_nn(:,target) = classes_obs_target_opt(idx_opt_nn(:, target),target);
    z_target_opt_nn(:,target) = z_cal(idx_opt_nn(:, target));
    pmf_diff_z_plus_z_nn(:,target) = pmf_z_target_given_neigh_opt(idx_opt_nn(:, target),target);     
end

weights_cont_or_obs_nn = zeros(n_neighbor_andor,length(z_cal));
weights_cont_and_obs_nn = zeros(n_neighbor_andor,length(z_cal)); 
normalized_weight_or_cont_obs_nn = zeros(n_neighbor_andor,length(z_cal));

for target = 1:length(z_cal)
    for obs = 1:n_neighbor_andor
        weights_cont_or_obs_nn(obs,target) =  best_w_OR(classes_obs_nn(obs,target)) + w_OR_slope(classes_obs_nn(obs,target)) * ( mat_euc_distance_xy_nn(obs,target) - edges_distance_classes_range(classes_obs_nn(obs,target)) ); %continuous OR weights
        weights_cont_and_obs_nn(obs,target) =  best_w_AND(classes_obs_nn(obs,target)) + w_AND_slope(classes_obs_nn(obs,target)) * ( mat_euc_distance_xy_nn(obs,target) - edges_distance_classes_range(classes_obs_nn(obs,target)) ); %continuous OR weights

    end
end

for target = 1:length(z_cal)
    for obs = 1:n_neighbor_andor
        normalized_weight_or_cont_obs_nn(obs,target) = weights_cont_or_obs_nn(obs,target) ./ sum(weights_cont_or_obs_nn(:,target)); %continuous OR weights
     end
end

% calculate AND_OR PMF for target 
w_alpha_beta = NaN(length(alpha)*length(beta),2);
i=1; %counter
for alpha_ = alpha  %for each set of alpha weights 
    for beta_ = beta %for each set of alpha weights 
        for target = 1 : length(z_cal) %for each target
            idx = [1:n_neighbor_andor]; %it ignores when i=j
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
probab_w_alpha_beta_obs = NaN(size(w_alpha_beta,1),length(z_cal)); %probability of z_obs (columns) in the correspondent set of AND_OR weights (row)
DKL_w_alpha_beta = NaN(size(w_alpha_beta,1),1); %DKL between the true Z and the predicted one. Each row corresponds to a set of weights
probab_AND_true_obs = ones(size(w_alpha_beta,1),length(z_cal));

for w = 1 : size(w_alpha_beta,1) %for each set of AND_OR weights
    prob_ = pmf_alpha_beta_nn(w,:);
    PMF_true = ones(1,size(z_cal,1));
    DKL_w_alpha_beta(w) = f_performance_prob(z_cal', prob_, PMF_true, edges_z);
end


idx_best_alpha_beta = min(find(DKL_w_alpha_beta==min(DKL_w_alpha_beta))); %em caso de empate, vamos pegar o que tem mais OR
DKL_score_mean_cal = DKL_w_alpha_beta(idx_best_alpha_beta);
best_alpha = w_alpha_beta(idx_best_alpha_beta,1);
best_beta = w_alpha_beta(idx_best_alpha_beta,2);
pmf_z_obs_AND_OR = pmf_alpha_beta_nn(idx_best_alpha_beta,:); %the best z PMFs 

% if find(DKL_w_alpha_beta==min(DKL_w_alpha_beta)) == 1
%      disp('The best AND_OR combination is a pure AND.'); % TODO! why is this bad?
% end

% step 8.1: Plot the performance of weights AND_OR TODO check xticks %TODO
% this plot!!!!
%analyzing the results
% xticks_ = cell(size(alpha));
% for i = 1:length(DKL_w_alpha_beta)
%     xticks_{i} = [ num2str(alpha(i)) '/' num2str(beta(i))];
% end
bin_centers_edges_z = [];
bin_centers_edges_z = edges_z(1:length(edges_z)-1) + binwidth_z/2;

clear pmf_alpha_beta_nn

%% Geo3
%step 2: distance of the target to its neighbors
% calculate the euclidean distance between target/grid (column) and
% their known neighbors/observation (row)
x_target = x(:);
y_target = y(:);
% z_target = z(:);

mat_euc_distance_obs_target_pred = NaN(length(x_cal),length(x(:))); % create empty vector ?????????????????????????????????

for target = 1:numel(x_target)
    for obs = 1:length(x_cal)
        mat_euc_distance_obs_target_pred(obs,target) = f_euclidean_dist(x_target(target), y_target(target), x_cal(obs), y_cal(obs));
    end

end

% step 3: classify the target neighbors
%according to the distance to the target, indentify the class of each neighbor 

classes_obs_target_pred = zeros(length(x_cal),length(x_target));

for target = 1:numel(x_target)
    for obs = 1:length(x_cal)
        for class_ = 1 : n_classes_range
            if mat_euc_distance_obs_target_pred(obs,target) > edges_distance_classes_range(class_) & mat_euc_distance_obs_target_pred(obs,target) <= edges_distance_classes_range(class_+1); % in case it is within the current lag class
                classes_obs_target_pred(obs,target) = class_; 
            end          
        end
%         if mat_euc_distance_obs_target_pred(obs,target) > edges_distance_classes_range(n_classes_range) %the last class take all distances 
%             classes_obs_target_pred(obs,target) = n_classes_range; %save the identified class
%         elseif mat_euc_distance_obs_target_pred(obs,target) == 0 %if the distance between points is zero, it will receive the first class identifier 
%             classes_obs_target_pred(obs,target) = 1; %
%         end    
        if mat_euc_distance_obs_target_pred(obs,target) == 0 %if the distance between points is zero, it will receive the first class identifier
           classes_obs_target_pred(obs,target) = 1; %
       end
     
    end    
end  

% step 3.1: Distribution of the Neighbors | Neighbor Cardinality | Connectivity distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%REVER
%plot the # of pairs of the GRID (target) by class 
n_pairs_by_class_pred = [];
for i = 1:n_classes_range
    n_pairs_by_class_pred(i) = sum(classes_obs_target_pred == i, 'all');
end

%plot the mean of the # of measurement neighbors by class
nmean_obs_by_class_pred = NaN(n_classes_range,1);
[n_obs_by_class_pred, edgge_class] = hist(classes_obs_target_pred(:,:), -0.5:1:n_classes_range); %getting the # of observations by class for each target (row)
for class_ = 1 : n_classes_range+1 %for each class
    nmean_obs_by_class_pred(class_) = mean(n_obs_by_class_pred(class_,:)); %calculate the mean of neighbors
end

%plot the neighborless observations
target_idx_zero_neigh_pred = [];
for target = 1 : length(x_target) %for each target
    if n_obs_by_class_pred(end,target) == sum(n_obs_by_class_pred(2:end,target))
        target_idx_zero_neigh_pred = [target_idx_zero_neigh_pred target];
    end
end

% finding the closest ones

% step 3.1 Distribution of the Neighbors | Neighbor Cardinality | Connectivity distribution
idx_nn = NaN(n_neighbor,length(z_cal));
mat_euc_distance_obs_target_pred_nn = NaN(n_neighbor,length(z_cal));
classes_obs_target_pred_nn = NaN(n_neighbor,length(z_cal));
z_cal_nn = NaN(n_neighbor,length(z_cal));

for target = 1 : length(x_target) %for each target identidy the closest neighbors
    [~,idx_] = sort(mat_euc_distance_obs_target_pred(:,target));
    idx_nn(:, target) = idx_(1:n_neighbor);
    mat_euc_distance_obs_target_pred_nn(:,target) = mat_euc_distance_obs_target_pred(idx_nn(:, target),target);
    classes_obs_target_pred_nn(:,target) = classes_obs_target_pred(idx_nn(:, target),target);
    z_cal_nn(:,target) = z_cal(idx_nn(:, target));
end


% step 4: Normalized weights
weights_or_pred_nn = zeros(n_neighbor,length(x_target));
normalized_weight_or_pred_nn = zeros(n_neighbor,length(x_target));
weights_or_cont_pred_nn = zeros(n_neighbor,length(x_target));
normalized_weights_or_cont_pred_nn = zeros(n_neighbor,length(x_target));
weights_and_pred_nn = zeros(n_neighbor,length(x_target));
weights_and_cont_pred_nn = zeros(n_neighbor,length(x_target));

for target = 1:length(x_target)
    for obs = 1:n_neighbor
        weights_or_pred_nn(obs,target) = best_w_OR(classes_obs_target_pred_nn(obs,target)); %discrete OR weights
        weights_or_cont_pred_nn(obs,target) =  best_w_OR(classes_obs_target_pred_nn(obs,target)) + w_OR_slope(classes_obs_target_pred_nn(obs,target)) * ( mat_euc_distance_obs_target_pred_nn(obs,target) - edges_distance_classes_range(classes_obs_target_pred_nn(obs,target)) ); %continuous OR weights
        weights_and_pred_nn(obs,target) = best_w_AND(classes_obs_target_pred_nn(obs,target)); %discrete OR weights
        weights_and_cont_pred_nn(obs,target) =  best_w_AND(classes_obs_target_pred_nn(obs,target)) + w_AND_slope(classes_obs_target_pred_nn(obs,target)) * ( mat_euc_distance_obs_target_pred_nn(obs,target) - edges_distance_classes_range(classes_obs_target_pred_nn(obs,target)) ); %continuous OR weights
    end
end

for target = 1:length(x_target)
    for obs = 1:n_neighbor
        normalized_weight_or_pred_nn(obs,target) = weights_or_pred_nn(obs,target) ./ sum(weights_or_pred_nn(:,target)); %discrete OR weights
        normalized_weights_or_cont_pred_nn(obs,target) = weights_or_cont_pred_nn(obs,target) ./ sum(weights_or_cont_pred_nn(:,target)); %continuous OR weights
     end
end

% steps 5 and 6: z PMF of the known observations(obs_diff_z + z) + z PMF with OR_AND combination
pmf_alpha_beta_pred_nn = cell(1,numel(x_target));
target_idx_without_AND_pred = []; %vector for the observations which do not support AND combination
i= 0;

%f_ = waitbar(0,'Obtaining z PMFs...');

pmf_z_target_given_neigh_pred_ = nan(n_neighbor,n_bins_z); %z PMF for each neighbor (row)
diff_n_bins_shift = n_bins_z*n_bins_z_per_bin_shift - n_bins_shift; %how many small bins we will need to add (before and/or after) to have "numbins_z_per_bin_shift" bins for each bin in the final z PMF

for target = 1:numel(x_target)%2463%4873
    %waitbar(target/length(x_target),f_,sprintf('Obtaining z PMFs... (%.2f%%)', target/length(x_target)*100));
    %disp([num2str(target/length(x_all)*100),' %'])
    % step 5: z PMF of the known observations(obs_diff_z + z)
    for obs = 1:n_neighbor
        class_ = classes_obs_target_pred_nn(obs,target); 

        bins_shift = round((edges_diff_z_shift(1) + z_cal_nn(obs,target) - edges_z(1)) / binwidth_shift);
        idx_bins_shift = reshape( ...           % index matrix to find which bins from the shift PMF go into each bin of the z PMF
              [ones(1,bins_shift) 2:n_bins_shift+1 ones(1,int16(diff_n_bins_shift) - bins_shift)], ... %fill edges with 1 (zero bin)
              int16(n_bins_z_per_bin_shift), n_bins_z); % matrix with one column per z PMF bin and "numbins_z_per_bin_shift" rows
%             [ones(1,bins_shift) 2:numbins_shift+1 ones(1,diff_numbins_shift - bins_shift)], ... %fill edges with 1 (zero bin)
%             numbins_z_per_bin_shift, numbins_z); % matrix with one column per z PMF bin and "numbins_z_per_bin_shift" rows
        pmf_ = pmf_diff_z_by_class_obs_range_shift(class_, :);
        pmf_z_target_given_neigh_pred_(obs,:) = sum(pmf_(idx_bins_shift)); %group small bins from "pmf_diff_z_by_class_obs_range_shift" into large z PMF bins
        pmf_z_target_given_neigh_pred_(obs,:) = pmf_z_target_given_neigh_pred_(obs,:) / sum(pmf_z_target_given_neigh_pred_(obs,:)); %normalize
        %         [pmf_diff_z_plus_z_pred_(obs,:), ~] = histcounts(zm(obs) + obs_diff_z_by_class_range{1,class}, edges_z,  'Normalization', 'probability');
    end

    % step 6: z PMF with OR_AND combination
    %predicting the weighted PMF_AND_OR for each xi,yi_analized    
    idx = [1:n_neighbor];
    pmfs_ = pmf_z_target_given_neigh_pred_(idx,:);
    weights_or_ = normalized_weights_or_cont_pred_nn(idx,target);
    weights_and_ = weights_and_cont_pred_nn(idx,target);
    pmf_and_ = f_loglinear_aggregation(pmfs_,weights_and_);
    pmf_or_ = f_linear_aggregation(pmfs_,weights_or_);
    pmf_ = f_loglinear_aggregation([pmf_and_;pmf_or_],[best_alpha;best_beta]);
%    pmf_ = f_mixpdfs_expAND_loglinear(pmfs_, best_alpha, best_beta, weights_or_, weights_and_); %0=OR & 1=AND
%     if (find(isnan(pmf_)))
%         i = i + 1;
%         pmf_ = f_mixpdfs(pmfs_, 0, weights_or_); %if we have a problem in the intersection AND, it will calculate just the OR combination
%         target_idx_without_AND_pred(i) = target; %save the index of the observation for analysis
%     end
    [ pmf_alpha_beta_pred_nn{1,target} ] = pmf_;
end

% plot uncertainty results
%calculate the entropy of the z PMFs for all predicted points 
H_z_pmf_by_class_pred = NaN(numel(x_target),1); %vector for the entropy of the grid (xi,yi_analized)
% H_z_pmf_by_class_obs = NaN(length(xm),1); %vector for the entropy of the known observations


for target = 1:numel(x_target) %for each target
    H_z_pmf_by_class_pred(target,1) = f_entropy(cell2mat(pmf_alpha_beta_pred_nn(1,target))); %calculate the entropy
end

z_target_entropy_pred_plot = reshape(H_z_pmf_by_class_pred, size(x,1), size(x,2));

% z mean, median, mode, probab
%predict the z with max probability in the z PMF 
z_target_mean_pred = NaN(1,numel(x_target));
z_target_median_pred = NaN(1,numel(x_target));
z_target_mode_pred = NaN(1,numel(x_target));

for target = 1:numel(x_target) %xi,yi_analized 
    %mean
    z_target_mean_pred(1,target) = (sum(bin_centers_edges_z .* (pmf_alpha_beta_pred_nn{1,target}))) / sum(pmf_alpha_beta_pred_nn{1,target}); %interpolated mean (not the bin center)
    
    %median
    cmf_median_ = cumsum(pmf_alpha_beta_pred_nn{1,target});
    idx_median_ = find(cmf_median_ >= 0.5, 1);
    x_median_ = [bin_centers_edges_z(idx_median_-1) bin_centers_edges_z(idx_median_) bin_centers_edges_z(idx_median_+1)];
    y_median_ = [cmf_median_(idx_median_-1) cmf_median_(idx_median_) cmf_median_(idx_median_+1)];
    z_target_median_pred(1,target) = interp1(y_median_,x_median_,0.5); %interpolated median (not the bin center)
    
    %mode
    [~,idx_] = max(cell2mat(pmf_alpha_beta_pred_nn(1,target)));
    z_target_mode_pred(1,target) = bin_centers_edges_z(idx_); %bin center
end

%plot
z_target_mean_pred_GRID_plot = reshape(z_target_mean_pred, size(x,1), size(x,2));
z_target_median_pred_GRID_plot = reshape(z_target_median_pred, size(x,1), size(x,2));
z_target_mode_pred_GRID_plot = reshape(z_target_mode_pred, size(x,1), size(x,2));

clims = [min([z(:);z_target_mean_pred_GRID_plot(:)]) max([z(:);z_target_mean_pred_GRID_plot(:)])];

figure;
subplot(2,2,1);
pcolor(x,y,z);
xlabel('x [m]');
ylabel('y [m]');
title('Random field | Original');
shading flat;
h=colorbar;
caxis(clims);

subplot(2,2,2);
pcolor(x,y,z_target_mean_pred_GRID_plot);
xlabel('x [m]');
ylabel('y [m]');
title('Random field | HER exp. AND (mean)');
shading flat;
h=colorbar;
caxis(clims);
% hold on;
% plot(xm,ym,'k+');
%buffer
% polyout = polybuffer([xm ym],'points',2000);
% plot(polyout);

subplot(2,2,3);
mesh(x,y,z);
xlabel('x [m]');
ylabel('y [m]');
title('Random field | Original');
shading flat;
h=colorbar;
caxis(clims);
zlim([min(z, [], 'all'), max(z, [], 'all')]);

subplot(2,2,4);
mesh(x,y,z_target_mean_pred_GRID_plot);
xlabel('x [m]');
ylabel('y [m]');
title('Random field | HER exp. AND (mean)');
shading flat;
% colormap spring(5);
h=colorbar;
caxis(clims);
zlim([min(z, [], 'all'), max(z, [], 'all')]);

% sgtitle(sprintf('Cal/Val HER exp. AND| %0.0f cal. obs. / %0.0f val. obs.', dim_cal, dim_val_6384));
% saveas(gcf,sprintf('randomfield_output_her_dimens_cal%0.0f_neigh%0.0f_plusProb_expAND.fig', dim_cal, n_neighbor));

% clear and save
clear x_median_ y_median_ cmf_median_ idx_median_ str str1 weights_ pmfs_ idx probab_z_obs_ prob_ w obs class center_right center_left i class obs bin_ weights_ pmfs_ idx probab_z_obs prob w obs class center_right center_left i class obs_pred pmf_diff_z_plus_z_pred_ f_ idx_

%% Calculate performance metrics
z_target_pred_GRID = z_target_mean_pred_GRID_plot;
z_target_pred_GRID(idx_cal)=NaN;
z_target_entropy_pred = z_target_entropy_pred_plot;
z_target_entropy_pred(idx_cal)=NaN;

%%% performance validation set
% Root mean square deviation RMSD (RMSE)
% Mean Error (ME)  
% Mean Absolute Error (MAE) L1 norm, robust parameter estimator
% Nash-Sutcliffe model efficiency (r2, coefficient of determination)
[error_sign_val, RMSE_val, ME_val, MAE_val, NSE_val] = f_performance_det(z_target_pred_GRID(idx_val),z(idx_val));

% scoring rule - DKL
PMF_true_val = ones(1,length(idx_val));
PMF_simulated_val = pmf_alpha_beta_pred_nn(idx_val);
% load('C:\Users\Stephanie\Google Drive\04.2 PhD\05.Geostatistics IT\01.Matlab\2.random field G.P. paper dataset\edges_z_bw0.2.mat');
% shift_edges = find(edges_z_bw0pt2==edges_z(1))-1;
% for i = 1:length(idx_val)
%     PMF_simulated_val{1,i} = [zeros(1,shift_edges),PMF_simulated_val{1,i}]; %just because I am using edges different from the PMF
% end
DKL_score_mean_val = f_performance_prob(z_val', PMF_simulated_val, PMF_true_val, edges_z);%edges_z_bw0pt2);
correl_val = corr(z_val, (z_target_pred_GRID(idx_val) - z_val));

%%% performance test set
% Root mean square deviation RMSD (RMSE)
% Mean Error (ME)  
% Mean Absolute Error (MAE) L1 norm, robust parameter estimator
% Nash-Sutcliffe model efficiency (r2, coefficient of determination)
[error_sign_test, RMSE_test, ME_test, MAE_test, NSE_test] = f_performance_det(z_target_pred_GRID(idx_test),z(idx_test));

% scoring rule - DKL
PMF_true_test = ones(1,length(idx_test));
PMF_simulated_test = pmf_alpha_beta_pred_nn(idx_test);
DKL_score_mean_test = f_performance_prob(z_test', PMF_simulated_test, PMF_true_test, edges_z);
correl_test = corr(z_test, (z_target_pred_GRID(idx_test) - z_test));
