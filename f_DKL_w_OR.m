function [DKL_w_OR, varargout] = f_DKL_w_OR(weights_or, z_target_opt, classes_obs, pmf_diff_z_plus_z, edges_z)
    % weights: set of weights for the classes
    % z_target_opt: true value of the target (used to calculate DKL between the prediction and the truth)
    % classes_obs: class of the neighbors (according to the target)
    % pmf_diff_z_plus_z:

    % step 5: Normalized weights
    %Normalization of the weights. For each target, they need to sum 1.
    weights_obs = zeros(size(classes_obs,1),size(classes_obs,2)); %weights of the neighbors, based on the lag class (neighbors, target)
    normalized_weight_obs = zeros(size(classes_obs,1),size(classes_obs,2)); %normalized weights of the neighbors, where each column sum to 1 (neighbors, target, set of weights)

    for target = 1 : size(classes_obs,2) %for each target
        for i = 1 : size(classes_obs,1) %for each neighbor
            if classes_obs(i,target) ~= 0 % in case there is an associated class
                weights_obs(i,target) = weights_or(classes_obs(i,target)); %save the weight of the class 
            end        
        end
    end

    for target = 1 : size(classes_obs,2) %for each target
        for i = 1 : size(classes_obs,1) %for each neighbor
            normalized_weight_obs(i,target) = weights_obs(i,target) ./ sum(weights_obs(:,target)); %normalize the weights of the column 
        end
    end

    % step 7: OR PMF - Cross-validation
    %predict the z PMF of the target based on the OR combination of the neighbors
    %cross-validation: obtain the probability of the z_obs on the predicted z PMF 
    %select the best set of weights for the OR combination based on the 

    pmf_OR = cell(1,size(classes_obs,2)); %cell for the predicted z PMFs 

    %calculate weighted z PMF_OR of the target
    for target = 1 : size(classes_obs,2) %for each target
        idx = [1:target-1 target+1:size(classes_obs,1)]; %it jumps when target = neighbor (when i=j)
        pmfs_ = cell2mat(pmf_diff_z_plus_z(idx,target)); %take the PMF
        weights_ = normalized_weight_obs(idx,target); %take the normalized weights
%         [ pmf_OR{1,target} ] = f_mixpdfs(pmfs_, 0, weights_); %save the predicted z PMF of the target (0 = purely OR combination)
        [ pmf_OR{1,target} ] = f_linear_aggregation(pmfs_, weights_); %save the predicted z PMF of the target (0 = purely OR combination)

    end

    PMF_true = ones(1,size(classes_obs,1));
    [DKL_w_OR] = f_performance_prob(z_target_opt, pmf_OR, PMF_true, edges_z);

    % DKL minimization: obtaind the DKL between the true value and the predicted PMF
%     probab_w_OR_obs = NaN(1,size(classes_obs,1)); %probability of z_obs (columns) in the correspondent set of OR weights (row)
%     probab_true_obs = ones(1,size(classes_obs,1));
%     DKL_OR_true_pred = NaN(1,size(classes_obs,1));
% 
%     for target = 1 : size(classes_obs,2) %for each target
%         for i = 1 : length(edges_z) %for each bin of z
%             if z_target_opt(1,target) >= edges_z(i) & z_target_opt(1,target) < edges_z(i+1) % in case z_target is within the current bin
%                 prob_ = cell2mat(pmf_OR(1,target)); %temporary vector of the PMF 
%                 probab_w_OR_obs(1,target) = prob_(1,i); %vector with the probability of z_target                 
%             end
%         end
%     end
% 
% 
%     for target = 1 : length(z_target_opt) %for each target
%         DKL_OR_true_pred(1,target) = (log2(probab_true_obs(1,target)) - log2(probab_w_OR_obs(1,target)))*probab_true_obs(1,target); %calculate the DKL between the true value and the prediction
%     end
% 
% 
%     DKL_w_OR(:,1) = mean(DKL_OR_true_pred, 2); %accumulates the DKL of all predictions

    if nargout >= 2
        varargout{1} = normalized_weight_obs;
        varargout{2} = pmf_OR;
    end
end
