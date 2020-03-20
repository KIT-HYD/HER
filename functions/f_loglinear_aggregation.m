function [ pdf_out ] = f_loglinear_aggregation(pdfs, wand_pdfs)
%% function for log-linear aggregation of distributions
% - Goal is to combine m pdf's to a single pdf.
% - AND combination: The result is the intersection of all pdf's
% - wand_pdfs allow a weighted combination of pdf's (exponent of input pdf's)

% -------------- Input --------------
% - pdfs       [m,n]   matrix of m pdfs to be combined. Each pdfs has n bins
% - wand_pdfs  [m,1]   vector with weights for each of the m pdf's

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
    if ~isempty(find(isnan(wand_pdfs)))
        error('NaNs in wand_pdfs.')
    end

    % check probabilities in 'pdfs' sum to 1
    if ~isempty(find(abs(sum(pdfs,2) - 1) > .00001)) 
        error('Probablities in pdfs dont sum to 1.')
    end

    % check for equal input dimensions
    if ~isequal(size(pdfs,1),size(wand_pdfs,1))
        error('# of rows in pdfs and wand_pdfs must be equal.')
    end

    % Case 1: AND-only --> we only need to create the pdf_AND
    [~, cols] = find(pdfs == 0); % find all cols (=bins) with at least one zero 
    dummy = pdfs;
    if sum(dummy(:)) > 0  
         pdf_AND = sum(log(dummy).* wand_pdfs,1); 
         pdf_AND = exp(pdf_AND - max(pdf_AND)); % avoid numerical error due to very small numbers
         pdf_AND = pdf_AND / sum(pdf_AND); % normalize to p-sum=1
    else  % if the pdfs have contradicted each other for each bin, 
        error('AND-combination resulted in empty set'),
    end

    pdf_out = pdf_AND; % assign the pdf_AND to the output pdf

end

