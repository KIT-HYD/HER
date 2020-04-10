function [error_sign, RMSE, ME, MAE, NSE] = f_performance_det(z_simulated, z_true)
    % z_simulated = vector of the simulated results
    % z_true = vector of the true values
    
    error_sign = (z_simulated - z_true);
    error_abs = abs(error_sign);
    error_2 = error_sign .^2;
    
    % Root mean square deviation RMSD
    RMSE = sqrt(sum(error_2) / length(error_2)); %root mean square deviation
    
    % Mean Error (ME)  
    ME = sum( error_sign ) / length(z_simulated);
    
    % Mean Absolute Error (MAE) L1 norm, robust parameter estimator
    MAE = sum( error_abs ) / length(z_simulated);
    
    % Nash-Sutcliffe model efficiency (r2, coefficient of determination)
    NSE = 1 - sum(error_2) / sum( (z_true - mean(z_true)).^2 );

end