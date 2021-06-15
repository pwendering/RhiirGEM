function [Eineq, bineq, Eeq, beq] = getEnzymeConstraints(model, kcats)
%% [Eineq, bineq, Eeq, beq] = getEnzymeConstraints(model, kcats)
% returns:
% - an inequality matrix of size (#rules + #AND-rules) x (#reactions + #genes)
% the enzyme constraints are set according to the GPR rules:
% case 1) A OR B:  v_i <= kcat_i * (E_g,A + E_g,B)
% case 2) A AND B: v_i <= kcat_i * (E_g,A)
%                  v_i <= kcat_i * (E_g,B)
% complex rules: split GPR rules by AND and add constraints each component
% Input:
%       struct model:       metabolic model
%       double kcats:       turnover values for every reaction
% Output:
%       double Eineq:       matrix with left hand side of inequality constraints
%       double bineq:       vector with right hand side of inequality constraints

% get numbers of genes and reactions
nGenes = numel(model.genes);
nRxns = size(model.S,2);

% small numbers alpha, beta, and epsilon for protein concentration limits
a = 1E-10; % lower bound for y==1

% initilize the reaction-protein pair counter
rxnProtPairCounter = 0;
% initialize reaction-protein association matrix
rxnProtMat = zeros(nRxns,nGenes);
% non-zero entries indicate an association of protein and reaction, the value
% of the entry indicated the index of the reaction-protein pair
for i=1:nRxns
    geneIdx = cellfun(@str2double,regexp(model.rules{i},'\d+','match'));
    rxnProtPairCounter = rxnProtPairCounter+numel(geneIdx);
    rxnProtMat(i,geneIdx) = rxnProtPairCounter-numel(geneIdx)+1:rxnProtPairCounter;
end

% ~~~~ Inequality matrix ~~~~ %

% initialize the inequality matrix (without AND-components)
andRules = cellfun(@(x)regexp(x, '&', 'match'), model.rules, 'un', 0);
andRules = [andRules{:}];
Eineq = zeros(sum(~cellfun(@isempty,model.rules))+numel(andRules),...
    nRxns + nGenes + rxnProtPairCounter + nGenes);
clear andRules

% initialize the matrix of constraints for for AND-connected genes, which are
% OR-connected (will increase in size)
andCompConstraints = zeros(1,size(Eineq,2));

% initialize row counter
rowIdx = 0;
andCompConstraintCount=0;
for i=1:nRxns
    gpr = model.rules{i};
    % only if the reaction has both a kcat value and a GPR rule assigned
    if kcats(i)>0 && ~isempty(gpr)
        
        if contains(gpr, '&')&&contains(gpr, '|')
            
            % if complex rule: split by main operator:
            
            % add white space on both sides of the rule
            gpr=[' ',gpr,' '];
            gpr = regexprep(gpr,'))',') )');
            gpr = regexprep(gpr,'(x', '( x');
            
            % complex rules have additional pairs of brackets, which is used to
            % split the rule into its major components; keep the matching separator
            [gprSplit,sep] = regexp(gpr,'[\|\&] \( ','split','match');
            splitChar = cellfun(@(x)regexp(x,'[\&\|]','match'),sep,'un',0);
            
            [gprSplit,sep] = regexp(gprSplit,' \) [\|\&]','split','match');
            sep = [sep{:}];
            splitChar = [splitChar,cellfun(@(x)regexp(x,'[\&\|]','match'),sep,'un',0)];
            
            % convert cell to cellstr
            gprSplit = [gprSplit{:}];
            splitChar = unique([splitChar{:}]);
            
            if numel(splitChar)>1
                error('assumption violated that only one main operator exists!')
            elseif numel(splitChar)==0
                error('rule recognized as complex but no main components have been found!')
            end
            
            if isequal(splitChar,{'|'})
                andCompConstraintCount = andCompConstraintCount+1;
                % if OR operator: add an additional variable for every
                % AND-connected component
                
                % initialize row with additional columns for complexes
                % (expands the matrix)
                row = zeros(1,size(andCompConstraints,2)+numel(gprSplit));
                
                % flux is lower than or equal to the sum of both
                % AND-connected components
                % v_i - kcat * (sum of alternative complexes) <= 0
                row(i) = 1;
                
                andCompIdx = size(andCompConstraints,2) + (1:numel(gprSplit));
                row(andCompIdx) = -kcats(i);
                
                andCompConstraints(andCompConstraintCount,1:max(andCompIdx)) = row;
                
                % concentration of AND-connected components is lower than
                % or equal to the minimum of the involved genes
                % E_gA&E_gB <= E_gA; E_gA&E_gB <= E_gB
                for j=1:numel(gprSplit)
                    % add a constraint for each of the partaking genes
                    % (reaction-specific)
                    EgIdx = rxnProtMat(i,cellfun(@str2double,regexp([gprSplit{j}], '\d+', 'match')));
                    EgIdx = nRxns + nGenes + EgIdx;
                    for k=1:numel(EgIdx)
                        andCompConstraintCount = andCompConstraintCount+1;
                        row = zeros(1,size(andCompConstraints,2));
                        % entry for the AND component
                        row(andCompIdx(j)) = 1;
                        % entry for the involved genes
                        row(EgIdx(k)) = -1;
                        % add row to andCompMatrix
                        andCompConstraints(andCompConstraintCount,:) = row;
                    end
                end

                % empty gprSplit so no additional constraints will be added
                % lateron
                gprSplit = [];
            else
                % if AND is main operator: upper bounds for each sum of OR
                % components
                gprSplit = strsplit(gpr,'&');
            end
            
        else
            gprSplit = strsplit(gpr,'&');
        end
        
        for j=1:numel(gprSplit)
            
            % ~~~~ v_i upper bound ~~~~ %
            % v_i - kcat * E_g(GPR_i) <= 0
            
            % initialize row
            rowIdx = rowIdx+1;
            row = zeros(1, size(Eineq,2));
            
            % find reaction-protein pair indices for the associated
            % proteins (current split)
            EgIdx = nRxns + nGenes + rxnProtMat(i,str2double(regexp(gprSplit{j}, '\d+', 'match')));
            % entry: sum of OR connected [E] * kcat
            row(EgIdx) = -kcats(i);
            % entry: v
            row(i) = 1;
            Eineq(rowIdx,:) = row;
            
        end
    end
end

% add enzyme promiscuity constraints [the sum of concentrations of an enzyme
% equals the total concentration of that protein across all reactions it is 
% catalyzing (or a protein complex it is participating in)]
% sum_i(E_g(i)) - E_g = 0
enzPromiscuityMatrix = zeros(numel(nGenes),size(Eineq,2));
for i=1:nGenes
    enzPromiscuityMatrix(i,nRxns+i) = -1;
    enzPromiscuityMatrix(i,nRxns+nGenes+rxnProtMat(any(rxnProtMat(:,i),2),i)) = 1;
end

% add boundaries for protein concentrations:
% ~~~~ E_g lower bounds ~~~~ %
% alpha*y - Eg <= 0
geneLowerBoundMatrix = [
    zeros(nGenes,nRxns),...                 v
    -eye(nGenes),...                        E_g
    zeros(nGenes,rxnProtPairCounter),...    RPP
    a*eye(nGenes)...                        y
];

% ~~~~ E_g upper bounds ~~~~ %
% E_g - y <= 0
geneUpperBoundMatrix = [
    zeros(nGenes,nRxns),...                 v
    eye(nGenes),...                         E_g
    zeros(nGenes,rxnProtPairCounter),...    RPP
    -eye(nGenes)...                  y
];

% append E_g constraints to inequalities matrix
Eineq = [Eineq;...
    geneLowerBoundMatrix;...                E_g lb
    geneUpperBoundMatrix...                 E_g ub
];                  
% construct bineq vector
bineq = [zeros(size(Eineq,1),1);...
    zeros(size(geneLowerBoundMatrix,1),1);...      E_g lb
    zeros(size(geneUpperBoundMatrix,1),1)...       E_g ub
];

% append AND-component constraints to inequality matrix
andCompNumber = size(andCompConstraints,2)-size(Eineq,2);
Eineq = [Eineq, zeros(size(Eineq,1),andCompNumber);
         andCompConstraints];

% construct combined right hand side vector for inequality constraints
bineq = [bineq; zeros(andCompConstraintCount,1)];

% construct combined equality matrix
Eeq = [enzPromiscuityMatrix, zeros(nGenes,andCompNumber)];

% define right hand side vector for equality constraints
beq = zeros(size(enzPromiscuityMatrix,1),1);

% remove all-zero rows from the matrix Eineq
nzRows = any(Eineq,2);
Eineq = Eineq(nzRows,:);
bineq = bineq(nzRows);

end