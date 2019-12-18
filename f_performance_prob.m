function [DKL_score_mean] = f_performance_prob(z_true, PMF_simulated, PMF_true, edges_z)
    % z_true = vector of the true values [1,n]
    % PMF_simulated = cell of the simulated PMFs {1,n}
    % PMF_true = cell of the true PMFs {1,n} OR vector of ones [1,n]
    % edges_z = z edges of the PMF

    if ~iscell(PMF_true) %true PMF = 1
        probab_obs_val_pred = NaN(1,length(z_true)); %probability of z_obs (columns)
        DKL_true_val_pred = NaN(1,length(z_true));

        for target = 1 : length(z_true) %for each target
            for i = 1 : length(edges_z) %for each bin of z
                if z_true(1,target) >= edges_z(i) & z_true(1,target) < edges_z(i+1) % in case z_target is within the current bin
                    prob_ = cell2mat(PMF_simulated(1,target)); %temporary vector of the PMF 
                    probab_obs_val_pred(1,target) = prob_(1,i);               
                end
            end
        end

        for target = 1 : length(z_true) %for each target
            DKL_true_val_pred(1,target) = (log2(PMF_true(1,target)) - log2(probab_obs_val_pred(1,target)))*PMF_true(1,target); %calculate the DKL between the true value and the prediction
        end

        DKL_score = sum(DKL_true_val_pred); %accumulates the DKL of all predictions
        DKL_score_mean = DKL_score/length(DKL_true_val_pred);
        
    elseif iscell(PMF_true) %TODO verify if it is working
        probab_obs_val_pred = NaN(1,length(z_true)); %probability of z_obs (columns)
        probab_obs_val_true = NaN(1,length(z_true)); %probability of z_obs (columns)
        DKL_true_val_pred = NaN(1,length(z_true));

        for target = 1 : length(z_true) %for each target
            for i = 1 : length(edges_z) %for each bin of z
                if z_true(1,target) >= edges_z(i) & z_true(1,target) < edges_z(i+1) % in case z_target is within the current bin
                    prob_ = cell2mat(PMF_simulated(1,target)); %temporary vector of the PMF 
                    probab_obs_val_pred(1,target) = prob_(1,i);               
                end
            end
        end

        for target = 1 : length(z_true) %for each target
            for i = 1 : length(edges_z) %for each bin of z
                if z_true(1,target) >= edges_z(i) & z_true(1,target) < edges_z(i+1) % in case z_target is within the current bin
                    prob_ = cell2mat(PMF_true(1,target)); %temporary vector of the PMF 
                    probab_obs_val_true(1,target) = prob_(1,i);               
                end
            end
        end
        
        for target = 1 : length(z_true) %for each target
            DKL_true_val_pred(1,target) = (log2(probab_obs_val_true(1,target)) - log2(probab_obs_val_pred(1,target)))*probab_obs_val_true(1,target); %calculate the DKL between the true value and the prediction
        end

        DKL_score = sum(DKL_true_val_pred); %accumulates the DKL of all predictions
        DKL_score_mean = DKL_score/length(DKL_true_val_pred);
    end
    
        
end