function [DKL_w_AND, varargout] = f_DKL_w_AND(weights_and, z_target_opt, classes_obs, pmf_diff_z_plus_z, edges_z)
    % weights: set of weights for the classes
    % z_target_opt: true value of the target (used to calculate DKL between the prediction and the truth)
    % classes_obs: class of the neighbors (according to the target)
    % pmf_diff_z_plus_z:

    % step 5: Normalized weights
    %Normalization of the weights. For each target, they need to sum 1.
%     weight_or_obs = zeros(size(classes_obs,1),size(classes_obs,2)); %weights of the neighbors, based on the lag class (neighbors, target)
%     normalized_weight_or_obs = zeros(size(classes_obs,1),size(classes_obs,2)); %normalized weights of the neighbors, where each column sum to 1 (neighbors, target, set of weights)
    weight_and_obs = zeros(size(classes_obs,1),size(classes_obs,2));
%     normalized_weight_and_obs = zeros(size(classes_obs,1),size(classes_obs,2));
    
    for target = 1 : size(classes_obs,2) %for each target
        for i = 1 : size(classes_obs,1) %for each neighbor
            if classes_obs(i,target) ~= 0 % in case there is an associated class
%                 weight_or_obs(i,target) = weights_or(classes_obs(i,target)); %save the weight of the class 
                weight_and_obs(i,target) = weights_and(classes_obs(i,target)); %save the weight of the class 
            end        
        end
    end

%     for target = 1 : size(classes_obs,2) %for each target
%         for i = 1 : size(classes_obs,1) %for each neighbor
%             normalized_weight_or_obs(i,target) = weight_or_obs(i,target) ./ sum(weight_or_obs(:,target)); %normalize the weights of the column 
%             %normalized_weight_and_obs(i,target) = weight_and_obs(i,target) ./ sum(weight_and_obs(:,target)); %normalize the weights of the column 
% 
%         end
%     end

    % step 7: OR PMF - Cross-validation
    %predict the z PMF of the target based on the OR combination of the neighbors
    %cross-validation: obtain the probability of the z_obs on the predicted z PMF 
    %select the best set of weights for the OR combination based on the 

    pmf_AND = cell(1,size(classes_obs,2)); %cell for the predicted z PMFs 

    %calculate weighted z PMF_OR of the target
    for target = 1 : size(classes_obs,2) %for each target
        idx = [1:target-1 target+1:size(classes_obs,1)]; %it jumps when target = neighbor (when i=j)
        pmfs_ = cell2mat(pmf_diff_z_plus_z(idx,target)); %take the PMF
%         weights_or_ = normalized_weight_or_obs(idx,target); %take the normalized weights
        weights_and_ = weight_and_obs(idx,target); %take the weights
%         [ pmf_AND{1,target} ] = f_mixpdfs_expAND(pmfs_, 1, weights_or_, weights_and_); %save the predicted z PMF of the target (0 = purely OR combination)
        [ pmf_AND{1,target} ] = f_loglinear_aggregation(pmfs_, weights_and_); %save the predicted z PMF of the target (0 = purely OR combination)
    end
    
    PMF_true = ones(1,size(classes_obs,1));
    [DKL_w_AND] = f_performance_prob(z_target_opt, pmf_AND, PMF_true, edges_z);

%     % DKL minimization: obtaind the DKL between the true value and the predicted PMF
%     probab_w_AND_obs = NaN(1,size(classes_obs,1)); %probability of z_obs (columns) in the correspondent set of OR weights (row)
%     probab_true_obs = ones(1,size(classes_obs,1));
%     DKL_AND_true_pred = NaN(1,size(classes_obs,1));
% 
%     for target = 1 : size(classes_obs,2) %for each target
%         for i = 1 : length(edges_z) %for each bin of z
%             if z_target_opt(1,target) >= edges_z(i) & z_target_opt(1,target) < edges_z(i+1) % in case z_target is within the current bin
%                 prob_ = cell2mat(pmf_AND(1,target)); %temporary vector of the PMF 
%                 probab_w_AND_obs(1,target) = prob_(1,i); %vector with the probability of z_target                 
%             end
%         end
%     end
% 
% 
%     for target = 1 : length(z_target_opt) %for each target
%         DKL_AND_true_pred(1,target) = (log2(probab_true_obs(1,target)) - log2(probab_w_AND_obs(1,target)))*probab_true_obs(1,target); %calculate the DKL between the true value and the prediction
%     end
% 
% 
%     DKL_w_AND(:,1) = mean(DKL_AND_true_pred, 2); %accumulates the DKL of all predictions
% 
    if nargout >= 2
        %varargout{1} = normalized_weight_or_obs;
        varargout{1} = pmf_AND;
    end
end
