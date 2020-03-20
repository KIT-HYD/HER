function [H] = f_entropy(pdf)
% Computes the Shannon Information Entropy of a discrete distribution in bit (log base 2)
% Input
% - pdf: [1,n] or [n,1] discrete probability density function (bin occupation frequencies, normalized to sum = 1)
%   Note
%   - pdf has to sum up to 1
%   - if pdf contains NaNs, the output will be H = NaN
% Output
% - H: Shannon entropy in [bit]
% Version
% - 2017/10/24 Uwe Ehret: handle the case of NaN's in the input
% - 2016/06/24 Uwe Ehret: intial version

% check if there are NaNs in pdf
if ~isempty(find(isnan(pdf)))
    H = NaN;
    return;
end

% check if pdf sums to 1
if abs(sum(pdf) - 1) > .00001
    error('Probablities dont sum to 1.')
end

% everything ok with the input --> compute H

    % initialize the output
    H = 0;

    % loop over all bins
    for i = 1 : length(pdf)
        if pdf(i) == 0 
            H = H;
        else
            H = H + (-log2(pdf(i))* pdf(i));   
        end
    end

end

