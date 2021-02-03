function [f, MW] = getAccountedProteinModel(model, uniProtFile, transcriptomicFile)
%% [f, mw] = getAccountedProteinModel(model, uniProtFile, transcriptomicFile)
% estimates the mass fraction of the products of the genes in a metabolic model
% Input:
%       struct model:               metabolic model
%       char uniProtFile:           path to UniProt file (at least with columns:
%                                   GeneNames, Sequence, (optional: Mass);
%                                   headers must match)
%       char transcriptomicFile:    path to file containing expression data
%                                   (tab or comma-separated; no header);
%                                   IDs in the first column should be
%                                   UniProt IDs or gene IDs
% Output:
%       double f:                   mass fraction
%       double MW:                  molecular weights for all proteins
%                                   contained in the model


% check number of arguments
if nargin < 3
    error('USAGE: [f, MW] = getAccountedProteinModel(model, uniProtFile, transcriptomicFile)')
end

% check input roughly
if ~isstruct(model) || ~isfield(model, 'genes')
    error('The model structure must have a field ''genes''.')
elseif ~ischar(uniProtFile) || ~exist(uniProtFile, 'file')
    error('UniProt file not found of input not of type char')
elseif ~ischar(transcriptomicFile) || ~exist(transcriptomicFile, 'file') || isempty(transcriptomicFile)
    disp('transcriptomic file not found of input not of type char')
    disp('')
    calcFractionBool=0;
else
    calcFractionBool=1;
end

tic
disp('Calculating the accounted mass fraction')
% amino acid molecular weights [g/mol]
MW_aa = [75.07,...  Glycine
    115.13,...      L-Proline
    89.09,...       L-Alanine
    117.15,...      L-Valine
    131.17,...      L-Leucine
    131.17,...      L-Isoleucine
    149.21,...      L-Methionine
    121.16,...      L-Cysteine
    165.19,...      L-Phenylalanine
    181.19,...      L-Tyrosine
    204.22,...      L-Tryptophan
    155.15,...      L-Histidine
    146.19,...      L-Lysine
    174.2,...       L-Arginine
    146.14,...      L-Glutamine
    132.12,...      L-Asparagine
    147.13,...      L-Glutamate
    133.1,...       L-Aspartate
    105.09,...      L-Serine
    119.12...       L-Threonine
    ];

% read UniProt file and select unique entries (wrt gene IDs)
uniprot = readtable(uniProtFile,'FileType','text',...
    'ReadVariableNames',true);
[geneIDs,ia]=unique(uniprot.GeneNames,'stable');

if ismember('Mass',uniprot.Properties.VariableNames)
    fprintf('\t> Using molecular weights from UniProt\n')
    MW = cellfun(@(x)str2double(regexprep(x,',','')),uniprot.Mass(ia));
    clear uniprot
else
    sequences = uniprot.Sequence(ia);
    clear uniprot
    
    %% determine amino acid compositions of all sequences
    fprintf('\t> Determining amino acid compositions for all proteins\n')
    aa_mat = zeros(numel(sequences), 20);
    for i=1:numel(sequences)
        coeff = getAminoAcidComposition(sequences{i});
        aa_mat(i,:) = coeff';
    end
    clear sequences coeff
    
    %% calculate total mass of all proteins
    disp('Calculating molecular weights')
    % multiply the occurrence of every amino acid with its respective molecular
    % weight
    MW = aa_mat * MW_aa';
    % correct for water
    MW = MW - 18.015*(sum(aa_mat,2)-1);
end

if calcFractionBool
    % read transcriptomic file
    fprintf('\t> Processing counts\n')
    transcriptomicTab = readtable(transcriptomicFile, 'ReadVariableNames', false);
    
    % find gene IDs in the UniProt file which are found in the transcriptomic file
    up2expressionMatch = ismember(geneIDs, transcriptomicTab.Var1);
    % take the subset
    geneIDs = geneIDs(up2expressionMatch);
    
    % find indices of gene ids in the transcriptomic file which are contained in
    % the UniProt file
    expression2upMatch = cellfun(@(x)find(ismember(transcriptomicTab.Var1, x)), geneIDs, 'un',0);
    expression2upMatch = cell2mat(expression2upMatch(~cellfun(@isempty,expression2upMatch)));
    
    % sort and subset the counts according to the IDs in the UniProt
    % file
    transcriptomicTab = transcriptomicTab(expression2upMatch,:);
    
    % scale relative to the sum
    transcriptomicTab.Var2 = transcriptomicTab.Var2 / sum(transcriptomicTab.Var2);
    
    % multiply the molecular weight with the protein count
    MW = MW(up2expressionMatch);
    
    ptotAll = MW' * transcriptomicTab.Var2;
    
    % calculate accounted mass fraction
    fprintf('\t> Calculating mass fraction\n')
    
    % find protein ids associated to the gene ids in the model
    up2modelMatch = cell2mat(cellfun(@(x)find(ismember(geneIDs,x)), model.genes, 'un', 0));
    
    fprintf('\t> %d genes could be associated to proteins in the UniProt file.\n', numel(up2modelMatch))
    
    % calculate the mass accounted by the model as above for the total mass
    ptotModel = MW(up2modelMatch)' * transcriptomicTab.Var2(up2modelMatch);
    
    % calculate the fraction
    f = ptotModel / ptotAll;
else
    up2modelMatch = cell2mat(cellfun(@(x)find(ismember(geneIDs,x)), model.genes, 'un', 0));
    f = 0.5;
end

% molecular weights
MW = MW(up2modelMatch);
fprintf('Finished! Time: %.2f min\n\n', toc/60)
end
