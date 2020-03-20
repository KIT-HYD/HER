function [ pdf_out ] = f_linear_aggregation(pdfs, wor_pdfs)%(pdfs, w_andor, w_pdfs)
% combines discrete pdfs
% Idea
% - Goal is to combine m pdf's to a single pdf. Possible combinations are 
%   AND: The result is the intersect of all pdfs
%   OR: The result is the combination of all pdfs
% - w_andor allows to mixture of the AND and OR combination
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
%   - OR: If we apply OR, we assume that there is some uncertainty in each
%     pdf-statement, but we assume that
%     - together they are honest about the value range (the OR combi will always yield pdfs
%       that at most cover the joint value range of all pdfs, but not more).
%     - All pdfs learn from the same underlying data set, but have not seen
%       all the data (OR essentially means we combine rows in a data set).
%       This also means that the precision of the resulting pdf is >= that
%       of each of the input pdfs
%   - When to use weights in the OR-combination?
%     This accounts for the different sizes of the samples underlying the
%     input pdf's. E.g one could have been constructed from 100
%     observations, the other from 1000. This needs to be considered by
%     giving the latter a weight of 1000/1100 and the first 100/1100.
%     Together the weights must sum to 1
%   - When to use weights between an AND and an OR combination?
%     If we are not fully certain that the different input pdfs are
%     completely honest, but still think together they no more than each of
%     them alone. So its a matter of consistency with the true underlying
%     system, and the effect of sample size and the dimensionality of the
%     system.
%   - How to expand value ranges for the resulting pdf?
%     For all of the above combinatons, the range of the resulting pdf will
%     never exceed that of all inputs, and for bins where all input pdfs
%     had p=0, the resulting pdf will also show p=0. If we think we suffer
%     from a limited data set, and states where p=0 in all inputs may be
%     observed later, we can add an uniform pdf spanning the entire
%     physically feasible value range, and with small but non-zero p's for
%     all bins. In an OR-combination, if the weigth of this uniform pdf is
%     non-zero, p will be non-zero in the output pdf.
%   - Why is it not possible to assign weights for the pdfs in a pure AND-combination?
% %     If we apply AND, we state that the different pdfs are completely
% %     consistent with the underlying truth, and that the different pdfs, as
% %     they have learned from the same data, can never contradict each
% %     other. In this case, if we see a weigth as a reflection of the # of
% %     data the pdf was constructed from, the weigths for all input pdfs
% %     must be the same (as they learned from the same # of data).
% - NOte: If all bins in all pdfs have non-zero p, and if weights for each pdf
%   are equal, then AND-only, OR-only and all possible combinations (values
%   of w_andor) yield the same result

% Input
% - pdfs: [m,n] matrix of m pdfs to be combined. Each pdfs has n bins
% - w_andor: [1,1]: Weight for the AND-combination. 
%                   Must be [0,1]. 0= no weight for AND-combination, full weight for OR-combination
% - w_pdfs: [m,1] vector with weights for each of the m pdf's. Must sum up to 1
% Output
% - pdf_out: [1,n] vector of probabilities of the mixed (merged) pdfs

% Version
% - 2018/07/03 Uwe Ehret: intial version

% check if there are NaNs in 'pdfs'
if ~isempty(find(isnan(pdfs)))
    error('NaNs in pdfs.')
end
% 
% % check if 'w_andor' is NaN
% if isnan(w_andor)
%     error('NaNs in w_andor')
% end

% check if there are NaNs in 'w_pdfs'
if ~isempty(find(isnan(wor_pdfs)))
    error('NaNs in wor_pdfs.')
end
% if ~isempty(find(isnan(wand_pdfs)))
%     error('NaNs in wand_pdfs.')
% end

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
% if ~isequal(size(pdfs,1),size(wand_pdfs,1))
%     error('# of rows in pdfs and wand_pdfs must be equal.')
% end

% % Case 1: AND-only --> we only need to create the pdf_AND
% if w_andor == 1
%     [~, cols] = find(pdfs == 0); % find all cols (=bins) with at least one zero 
%     %indx_badcols = unique(cols); % --> as at least one pdf said 'impossible', overall p for that bin will be zero
%     dummy = pdfs;
%     %dummy(:,indx_badcols) = 0;   % set the whole col to zeros %Fixed ST 
%     if sum(dummy(:)) > 0 %BEFORE: if ~isempty(dummy) %Fixed ST, since isempty does not catch zeros 
%         pdf_AND = sum(log(dummy).* wand_pdfs,1); %BEFORE: pdf_AND = sum(dummy,1);  % for the remaining non-zero cols, add the p's from all pdfs %Fixed ST, since AND should multiply
% %          pdf_AND = exp(pdf_AND); %ST
%          pdf_AND = exp(pdf_AND - max(pdf_AND)); % avoid numerical error due to very small numbers
%          pdf_AND = pdf_AND / sum(pdf_AND); % normalize to p-sum=1
%     else                         % if the pdfs have contradicted each other for each bin, something is wrong
%         error('AND-combination resulted in empty set'),
%     end
%     
%     pdf_out = pdf_AND; % assign the pdf_AND to the output pdf
    
% Case 2: OR-only --> we only need to create the pdf_OR
% elseif w_andor == 0
    dummy = pdfs .* wor_pdfs; % multiply each pdf with its weight
    pdf_OR = sum(dummy,1); % for each bin, sum up weighted p's from all pdfs
    pdf_OR = pdf_OR / sum(pdf_OR); % normalize to p-sum=1. This is actually not necessary
                                   % as we checked that the sum of
                                   % weights=1 and sum p for each pdf=1
                                   
    pdf_out = pdf_OR; % assign the pdf_OR to the output pdf                               
                                   
% Case 3: AND-OR combination --> we need to create both the pdf_AND and pdf_OR
% else
%     
%     % create the pdf_AND
%     [~, cols] = find(pdfs == 0); % find all cols (=bins) with at least one zero 
%     %indx_badcols = unique(cols); % --> as at least one pdf said 'impossible', overall p for that bin will be zero
%     dummy = pdfs;
%     %dummy(:,indx_badcols) = 0;   % set the whole col to zeros %Fixed ST 
%     if sum(dummy(:)) > 0 %BEFORE: if ~isempty(dummy) %Fixed ST, since isempty does not catch zeros 
%         pdf_AND = sum(log(dummy).* wand_pdfs,1); %BEFORE: pdf_AND = sum(dummy,1);  % for the remaining non-zero cols, add the p's from all pdfs %Fixed ST, since AND should multiply
% %         pdf_AND = exp(pdf_AND); %ST
%         pdf_AND = exp(pdf_AND - max(pdf_AND)); % avoid numerical error due to very small numbers
%         pdf_AND = pdf_AND / sum(pdf_AND); % normalize to p-sum=1
%     else                         % if the pdfs have contradicted each other for each bin, something is wrong
%         error('AND-combination resulted in empty set'),
%     end    
%     
%     % create the pdf_OR
%     dummy = pdfs .* wor_pdfs; % multiply each pdf with its weight
%     pdf_OR = sum(dummy,1); % for each bin, sum up weighted p's from all pdfs
%     pdf_OR = pdf_OR / sum(pdf_OR); % normalize to p-sum=1. This is actually not necessary
%                                    % as we checked that the sum of
%                                    % weights=1 and sum p for each pdf=1
%     
%     % combine pdf_AND and pdf_OR
%     pdf_out = (pdf_AND * w_andor) + (pdf_OR * (1-w_andor));
%     pdf_out = pdf_out / sum(pdf_out); % normalize to p-sum=1. This is actually not necessary
%     
% end
% 

end

