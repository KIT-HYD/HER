function [ pdf_out ] = f_loglinear_aggregation(pdfs, wand_pdfs)%(pdfs, w_andor, w_pdfs)
% combines discrete pdfs
% Idea
% - Goal is to combine m pdf's to a single pdf. Possible combinations are 
%   AND: The result is the intersect of all pdfs
% - w_pdfs allow, for the OR-combination, a weighted combination of pdfs
% - When to use which?
%   - AND: This means the results satisfies all pdfs at the same time. When
%     using AND, we assume that each pdf is honest (no uncertainty). This
%     implies that all these pdfs are consistent with the underlying true
%     data, but have maybe not absorbed all information from it. E.g one pdf
%     could state that the roll of a dice is even/odd, another could state
%     the roll is <3.5/>3.5. When we know both are honest, we can combine
%     them with AND to narrow down the range of possible values
%     Using AND means we assume each pdf has been learned from ALL the
%     data, but did not make exhaustive use of it. So its like adding
%     predictor columns in a data set, which means precision of the
%     resulting pdf is <= that of the each of the input pdfs (info can't
%     hurt)
%   - How to expand value ranges for the resulting pdf?
%     For all of the above combinatons, the range of the resulting pdf will
%     never exceed that of all inputs, and for bins where all input pdfs
%     had p=0, the resulting pdf will also show p=0. If we think we suffer
%     from a limited data set, and states where p=0 in all inputs may be
%     observed later, we can add an uniform pdf spanning the entire
%     physically feasible value range, and with small but non-zero p's for
%     all bins. In an OR-combination, if the weigth of this uniform pdf is
%     non-zero, p will be non-zero in the output pdf.

% Input
% - pdfs: [m,n] matrix of m pdfs to be combined. Each pdfs has n bins
% - w_pdfs: [m,1] vector with weights for each of the m pdf's. Must sum up to 1
% Output
% - pdf_out: [1,n] vector of probabilities of the mixed (merged) pdfs

% Version
% - 2018/07/03 Uwe Ehret: intial version

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
    %indx_badcols = unique(cols); % --> as at least one pdf said 'impossible', overall p for that bin will be zero
    dummy = pdfs;
    %dummy(:,indx_badcols) = 0;   % set the whole col to zeros %Fixed ST 
    if sum(dummy(:)) > 0 %BEFORE: if ~isempty(dummy) %Fixed ST, since isempty does not catch zeros 
        pdf_AND = sum(log(dummy).* wand_pdfs,1); %BEFORE: pdf_AND = sum(dummy,1);  % for the remaining non-zero cols, add the p's from all pdfs %Fixed ST, since AND should multiply
%          pdf_AND = exp(pdf_AND); %ST
         pdf_AND = exp(pdf_AND - max(pdf_AND)); % avoid numerical error due to very small numbers
         pdf_AND = pdf_AND / sum(pdf_AND); % normalize to p-sum=1
    else  % if the pdfs have contradicted each other for each bin, 
        error('AND-combination resulted in empty set'),
    end
    
    pdf_out = pdf_AND; % assign the pdf_AND to the output pdf
    
end

