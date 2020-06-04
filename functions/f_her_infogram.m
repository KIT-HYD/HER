function [her] = f_her_infogram(x_cal, y_cal,z_cal, her)
%% function to extract the Spatial structure of the data

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

% Geo1: Spatial Characterization
% step 1: get diff_z as a function of the distance 
% calculate the euclidean distance and the difference between all point pairs
zm_res_ = min(diff(unique(z_cal))); % smallest Z resolution
her.n_bins_z_per_bin_shift = min(100, round(0+her.binwidth_z / zm_res_)); % limit to 100 or binwidth_shift = 1 (in the case where all zm are integers)
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
her.edges_distance_classes = [0:her.lag_dist:her.lag_dist*her.n_lag]; %edges of the distance classes or DEFINE  it
% her.edges_distance_classes = [0,0.0250000000000000,0.0500000000000000,0.100000000000000,0.150000000000000,0.200000000000000,0.250000000000000,0.300000000000000,0.350000000000000,0.400000000000000,0.450000000000000,0.500000000000000,0.600000000000000,0.650000000000000,0.700000000000000,0.750000000000000,0.800000000000000,0.850000000000000,0.900000000000000,0.950000000000000,1,1.05000000000000,1.10000000000000,1.15000000000000,1.20000000000000,1.25000000000000,1.30000000000000,1.35000000000000,1.40000000000000,1.45000000000000,1.50000000000000,1.55000000000000,1.60000000000000,1.65000000000000,1.70000000000000,1.75000000000000,1.80000000000000,1.85000000000000,1.90000000000000,1.95000000000000,2,2.05000000000000,2.10000000000000,2.15000000000000,2.20000000000000,2.25000000000000,2.30000000000000,2.35000000000000,2.40000000000000,2.45000000000000,2.50000000000000,2.55000000000000,2.60000000000000,2.65000000000000,2.70000000000000,2.75000000000000,2.80000000000000,2.85000000000000,2.90000000000000,2.95000000000000,3,3.05000000000000,3.10000000000000,3.15000000000000,3.20000000000000,3.25000000000000,3.30000000000000,3.35000000000000,3.40000000000000,3.45000000000000,3.50000000000000,3.55000000000000,3.60000000000000,3.65000000000000,3.70000000000000,3.75000000000000,3.80000000000000,3.85000000000000,3.90000000000000,3.95000000000000,4,4.05000000000000,4.10000000000000,4.15000000000000,4.20000000000000,4.25000000000000,4.30000000000000,4.35000000000000,4.40000000000000,4.45000000000000,4.50000000000000,4.55000000000000,4.60000000000000,4.65000000000000,4.70000000000000,4.75000000000000,4.80000000000000,4.85000000000000,4.90000000000000,4.95000000000000,5,5.05000000000000,5.10000000000000,5.15000000000000,5.20000000000000,5.25000000000000,5.30000000000000,5.35000000000000,5.40000000000000,5.45000000000000,5.50000000000000,5.55000000000000,5.60000000000000,5.65000000000000];
% her.bin_centers_distance_classes = her.edges_distance_classes(1:end-1) + her.lag_dist/2; %UPDATED bin centers of the distance classes
her.bin_centers_distance_classes = NaN(1,her.n_lag);%bin centers of the distance classes
% check empty classes
count = 1;
her.emptyclass = [];
for i = 1 : her.n_lag %for each class
    idx = her.mat_euc_distance_xy > her.edges_distance_classes(i) & her.mat_euc_distance_xy <= her.edges_distance_classes(i+1); %find the diff_z within the current lag class
    her.bin_centers_distance_classes(i) = (her.edges_distance_classes(i) + her.edges_distance_classes(i+1))/2;
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
    her.n_classes_limit = 20;  %DEFINE it looking at the infogram
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

% For Geo2: define edges of Z 
max_z = max(z_cal(:));
min_z = min(z_cal(:));
max_diff_z = her.edges_diff_z_shift(end);
min_diff_z = her.edges_diff_z_shift(1);

% compute bin edges
n_bins_left = -floor((min_diff_z + min_z)/ her.binwidth_z);
n_bins_right = floor((max_diff_z + max_z)/ her.binwidth_z + 1); %edge <= x < edge
her.n_bins_z = n_bins_left+n_bins_right;
mini = -her.binwidth_z*n_bins_left;
maxi =  her.binwidth_z*n_bins_right;
her.edges_z = linspace(mini,maxi,her.n_bins_z+1); %edges of the z+diff_z pmf
[pmf_z_all_obs,~] = histcounts(z_cal, her.edges_z,'Normalization', 'probability'); %calculate the z PMF of the full dataset 

her.bin_centers_edges_z = [];
her.bin_centers_edges_z = her.edges_z(1:length(her.edges_z)-1) + her.binwidth_z/2;
her.bin_centers_edges_diff_z = [];
her.bin_centers_edges_diff_z = her.edges_diff_z(1:length(her.edges_diff_z)-1) + her.binwidth/2;
end