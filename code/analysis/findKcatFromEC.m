function kcatAll = findKcatFromEC(model, substrates, kcatFileName, taxTerm)
%% kcatAll = findKcatFromEC(model, substrates, kcatFileName, taxTerm)
% This function scans a file containing turnover numbers for matching 
% E.C. numbers (assigned to reactions) in a metabolic model. 
% The results from matching an E.C. number are, if possible, further 
% filtered for matching substrates or given taxonomic classification (any level).
% If kcat(s) have been found, the overall maximum is returned chosen per
% reaction (per E.C. number and over all associated E.C. numbers).
% Compounds, which are often non-specific are excluded from the matching
% procedure.
% Input:
%   struct model:           COBRA metabolic model structure
%   cell substrates:        array containing substrates per reaction
%   char kcatFileName:      path to file containing turnover numbers
%                           columns: E.C.-Number;substrate names;taxonomic
%                           lineage;turnover value;'*'
%                           (format as used in GECKO (Sanchez et al. 2017))
%                           the file has no header line!
%   char taxTerm:           taxonomic classification, can be organism name
%                           but also higher-level terms
% Output:
%   double kcatAll:         array containing one kcat per matched reaction;
%                           if unassigned: kcatAll(i)==0

% compounds which should be excluded from substrate matching
% exclude = {'nadh', 'nadp', 'nadph', 'nad', 'atp', 'adp', 'h2o', 'coa', 'h',...
%     'fadh2', 'fad'};
exclude={''};

% initialize outpur variable and counters
kcatAll = zeros(size(model.rxns));
ecNotFoundCount = 0;
prunedECMatchCount = 0;
substrateMatchCount = 0;
fullTaxMatchCount = 0;
fullMatchCount = 0;
% regex pattern to match E.C. numbers
ecPattern = '\d+\.\d+\.\d+\.\d+';

% determine the total number of E.C. numbers in the model
tmpEC = regexp(model.rxnECNumbers, ecPattern, 'match');
nEC  = sum(~cellfun(@isempty,tmpEC)); clear tmpEC

% loop over all reactions
n = numel(model.rxns);
fprintf('\n')
disp(['Obtaining maximum kcats for ', num2str(n), ' reactions.'])
fprintf('\n')
for i=1:n
    if mod(i,100)==0
        disp(['Processed ', num2str(i), ' reactions ...'])
    end
    
    % get all E.C. numbers assigned to the current reaction
    ec = regexp(model.rxnECNumbers{i}, ecPattern, 'match');
    
    % only proceed if the reaction has at least one E.C. number assigned
    if ~isempty(ec)
        % initialize kcat array for the current reaction
        kcat = zeros(size(ec));
        for j=1:numel(ec)
            % initial E.C. number level and default status
            level = 4;
            status = 1;
            % define the initial query
            query = ['^' ec{j} '\>'];
            
            while status ~= 0 && level ~= 1
                
                % call unix grep to scan turnover value file
                [status, res] = unix(['grep "', query, '" ', kcatFileName]);
                
                % check whether we have a zero exit code
                if status == 0 && ~isempty(res)
                    
                    % process the result to obtain a table
                    res = strtrim(res);
                    res = strsplit(res, {'\t', '\n'});
                    
                    % a check if every row was five column entries
                    try
                        res = cell2table(reshape(res, 5, numel(res)/5)');
                    catch
                        disp(res)
                    end
                    
                    % convert substrate names to a universal format
                    s_db = regexprep(res.(2), '^\(', '');
                    s_db = regexprep(s_db, '[\ \-\+\(\)\:\;\,\[\]]+', '_');
                    s_db = regexprep(s_db, '__', '_');
                    s_db = regexprep(s_db, '_$', '');
                    s_db = lower(regexprep(s_db, '^_', ''));
                    
                    % exclude unspecific compounds
                    substrates{i} = lower(substrates{i});
                    
                    for k=1:numel(exclude)
                        substrates{i} = regexprep(substrates{i}, ['^',exclude{k}, '$'], '');
                    end
                    
                    % find matching substrates
                    substrateMatchIdx=cellfun(@(x)any(ismember(strsplit(x,'|'),substrates{i})),s_db);
                    
                    % find entries for the given taxonomic term
                    taxMatchIdx = ~cellfun('isempty', regexpi(res.(3), taxTerm));
                    
                    if any(substrateMatchIdx&taxMatchIdx)
                        substrateMatchCount = substrateMatchCount + 1;
                        fullTaxMatchCount = fullTaxMatchCount + 1;
                        fullMatchCount = fullMatchCount + 1;
                        % substrate and fungal match
                        k = res.(4)(substrateMatchIdx&taxMatchIdx);
                    elseif any(substrateMatchIdx)
                        substrateMatchCount = substrateMatchCount + 1;
                        % only substrate match
                        k = res.(4)(substrateMatchIdx);
                    elseif any(taxMatchIdx)
                        fullTaxMatchCount = fullTaxMatchCount + 1;
                        % only taxonomic match
                        k = res.(4)(taxMatchIdx);
                    else
                        % no matches, take all
                        k = res.(4);
                    end
                    
                    if level < 4
                        prunedECMatchCount = prunedECMatchCount + 1;
                    end
                else
                    % remove the lower EC level
                    query = regexprep(query, '\.[\d\-]+(\\>)*$', '');
                    level = level - 1;
                end
                
            end
            
            if status ~= 0
                k = {'0'};
                if ~isempty(regexp(ec(j),ecPattern))
                    % only count as 'not found' if the reaction has an E.C.
                    % number assigned
                    ecNotFoundCount = ecNotFoundCount + 1;
                end
            end
            k = cellfun(@str2double, k);
            kcat(j) = max(k);
        end
        kcatAll(i) = max(kcat);
    end
    
end


fprintf('\nfinished!\n')
disp(['number of complete matches: ', num2str(fullMatchCount)])
disp(['total number of E.C. numbers: ', num2str(nEC)])
disp(['substrate matches: ', num2str(substrateMatchCount)])
disp(['taxonomic matches: ', num2str(fullTaxMatchCount)])
disp(['E.C. class matches: ', num2str(prunedECMatchCount)])
disp([num2str(ecNotFoundCount), ' E.C. numbers could not be matched'])
disp('')
end
