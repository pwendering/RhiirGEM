function res = splitGPRrule(gpr, val)
% splits GPR rule recursively, returns the sum for OR rules and
% minimum for AND rules
% Input:
%   char gpr:           GPR rule
%   double val:         values (e.g. expression) for all genes or proteins
% Output:
%   double res:         results of applying the rules for AND and OR
%                       operatators (minimum, sum) recursively

% get the gene indices
gene_idx = regexp(gpr, '\d+', 'match');
gene_idx = cellfun(@str2num, gene_idx);

% comlex rule containing AND and OR
if contains(gpr, '&')&&contains(gpr, '|')
    
    % split by main operator:
    
    % add white space on both sides of the rule
    gpr=[' ',gpr,' '];
    
    % complex rules have additional pairs of brackets, which is used to
    % split the rule into its major components; keep the matching separator
    [gprSplit,sep]=regexp(gpr,'[\|\&] \( ','split','match');
    splitChar=cellfun(@(x)regexp(x,'[\&\|]','match'),sep,'un',0);
    
    [gprSplit,sep]=regexp(gprSplit,' \) [\|\&]','split','match');
    sep=[sep{:}];
    splitChar=[splitChar,cellfun(@(x)regexp(x,'[\&\|]','match'),sep,'un',0)];
    
    % convert cell to cellstr
    gprSplit=[gprSplit{:}];
    splitChar=unique([splitChar{:}]);
    
    if numel(splitChar)>1
        error('assumption violated that only one main operator exists!')
    elseif numel(splitChar)==0
        error('rule recognized as complex but no main components have been found!')
    end
    
    if isequal(splitChar,{'|'})
        % if OR rule, take maximum
        res = sum(cellfun(@(x)splitGPRrule(x,val),gprSplit));
    else
        % if AND rule, take minimum
        res = min(cellfun(@(x)splitGPRrule(x,val),gprSplit));
    end
   
% rule contains only AND
elseif contains(gpr, '&')
    % take the minimum over all values
    res = min(val(gene_idx));
    
% singular or only OR
else
    res = sum(val(gene_idx));
end

end