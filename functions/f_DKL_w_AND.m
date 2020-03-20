function [DKL_w_AND, varargout] = f_DKL_w_AND(weights_and, z_target_opt, classes_obs, pmf_diff_z_plus_z, edges_z)
%% function for calculating DKL with AND aggregation of pdf's
% - Using leave-one-out cross-validation, the goal is to calculate the 
% mean DKL of the left-out targets using AND aggregation method.

% -------------- Input -------------- 
% - weights_and          [c,1]   vector with weights for each distance class. 
% - z_target_opt        [1,T]   true value of the z target (used to calculate 
%                               DKL between the prediction and the truth)
% - classes_obs         [T,T]   matrix with the class of observations (row)
%                               in relation to its distance to target (column). 
%                               Optimization: class zero when the neighbor = target. 
% - pmf_diff_z_plus_z   {T,T}   cell of T z_pdfs of the observations (rows) 
%                               to be combined for each target (column). 
%                               Each pdfs has n bins
% - edges_z             [1,n+1] edges of the bins of the z_PMF 

% -------------- Output --------------
% - DKL_w_AND                 {1,T}   mean DKL of the pdf aggregation
% - weight_and_obs            [T,T]   matriz with the weights contribution
%                                     of the observations to each target.
%                                     observation (row), target (column)
% - pmf_AND                   {1,T}   prediction of the pdf of each target
%                                     using the AND aggregation method

% -------------- Version --------------
% - 2019/10/01 Stephanie Thiesen: intial version
% - 2020/03/20 Stephanie Thiesen: input + output definition

% -------------- Script --------------
% Step 1: Associate the weights of observations according to its class in
% relation to the target
    weight_and_obs = zeros(size(classes_obs,1),size(classes_obs,2));
    for target = 1 : size(classes_obs,2) %for each target
        for i = 1 : size(classes_obs,1) %for each neighbor
            if classes_obs(i,target) ~= 0 % in case there is an associated class
                weight_and_obs(i,target) = weights_and(classes_obs(i,target)); %save the weight of the class 
            end        
        end
    end

% Step 2: AND PMF - Cross-validation
    %predict the z PMF of the target based on the AND combination of the neighbors
    %cross-validation: obtain the probability of the z_obs on the predicted z PMF 
    pmf_AND = cell(1,size(classes_obs,2)); %cell for the predicted z PMFs 
    %calculate weighted z PMF_AND of the target
    for target = 1 : size(classes_obs,2) %for each target
        idx = [1:target-1 target+1:size(classes_obs,1)]; %it jumps when target = neighbor (when i=j)
        pmfs_ = cell2mat(pmf_diff_z_plus_z(idx,target)); %take the PMF
        weights_and_ = weight_and_obs(idx,target); %take the weights
        [ pmf_AND{1,target} ] = f_loglinear_aggregation(pmfs_, weights_and_); %save the predicted z PMF of the target (0 = purely OR combination)
    end
    PMF_true = ones(1,size(classes_obs,1));
    [DKL_w_AND] = f_performance_prob(z_target_opt, pmf_AND, PMF_true, edges_z);

    if nargout >= 2
        varargout{1} = weight_and_obs;
        varargout{2} = pmf_AND;
    end
end
