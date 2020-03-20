function [ pdf_out ] = f_linear_aggregation(pdfs, wor_pdfs)
%% function for linear aggregation of distributions
% - Goal is to combine m pdf's to a single pdf.
% - OR combination: The result is the average of all pdf's (mixture distribution)
% - wand_pdfs allow a weighted combination of pdf's (exponent of input pdf's)

% -------------- Input --------------
% - pdfs       [m,n]   matrix of m pdfs to be combined. Each pdfs has n bins
% - wor_pdfs   [m,1]   vector with weights for each of the m pdf's. Must sum 1

% -------------- Output --------------
% - pdf_out [1,n]   vector of probabilities of the aggregated pdf's

% -------------- Version --------------
% - 2018/07/03 Uwe Ehret: intial version
% - 2020/03/20 Stephanie Thiesen: corrections and adjustments

% -------------- Script --------------
    % check if there are NaNs in 'pdfs'
    if ~isempty(find(isnan(pdfs)))
        error('NaNs in pdfs.')
    end

    % check if there are NaNs in 'w_pdfs'
    if ~isempty(find(isnan(wor_pdfs)))
        error('NaNs in wor_pdfs.')
    end

    % check probabilities in 'pdfs' sum to 1
    if ~isempty(find(abs(sum(pdfs,2) - 1) > .00001)) 
        error('Probablities in pdfs dont sum to 1.')
    end

    % check probabilities in 'wor_pdfs' sum to 1
    if abs(sum(wor_pdfs) - 1) > .00001
        error('Probablities in w_pdfs dont sum to 1.')
    end

    % check for equal input dimensions
    if ~isequal(size(pdfs,1),size(wor_pdfs,1))
        error('# of rows in pdfs and wor_pdfs must be equal.')
    end

    % Case 2: OR-only --> we only need to create the pdf_OR
    dummy = pdfs .* wor_pdfs; % multiply each pdf with its weight
    pdf_OR = sum(dummy,1); % for each bin, sum up weighted p's from all pdfs
    pdf_OR = pdf_OR / sum(pdf_OR); % normalize to p-sum=1. This is actually not necessary
                                   % as we checked that the sum of
                                   % weights=1 and sum p for each pdf=1

    pdf_out = pdf_OR; % assign the pdf_OR to the output pdf  
end

