function [growth, solution, ecModel] = simulateGrowthVmax(model, rxnKcats,...
    transcriptomicFile, C, f, MW)
%% modifies upper bound for reactions as kcat * [E] (factor 'f' from Sánchez et al. 2017)
% Input:
%   struct model:                   metabolic model structure
%   double rxnKcats:                turnover values for every reaction
%   char transcriptomicFile:        path to transcript counts per gene
%                                   (two-column csv/tab with columns GeneID,
%                                   Count without a header line)
%   double C:                       protein content [g/gDW]
%   double f:                       mass fraction which is accounted for
%                                   in the model
%   double MW:                      vector of molecular weights for each
%                                   gene product [g / mol]
% Output:
%   double growth:                  growth values at every concentration
%   double solution:               FBA results at different concentrations
%   struct ecModel:                   constrained model
if nargin<6
    error('Not enough input arguments')
end
ecModel = model;

%% process transcriptomic data
disp('Processing transcript counts')
% read file
transcriptCounts = readtable(transcriptomicFile);
% take the subset of genes in the model
idx = cell2mat(cellfun(@(x)find(ismember(transcriptCounts.Var1, x)), model.genes, 'un', 0));
transcriptCounts = transcriptCounts.Var2(idx);
% scale counts to get relative abundances
transcriptCounts = transcriptCounts / sum(transcriptCounts);
% assign the minimum value to zero entries to avoid complete knock-outs
transcriptCounts(transcriptCounts==0) = min(transcriptCounts(transcriptCounts>0));


% estimate protein concentrations from transcript abundance
transcriptCounts = ( 1000 * transcriptCounts * f * C ) ./  MW ;

% assign enzyme concentration to reactions according to the gpr rules
E = zeros(size(ecModel.S,2),1);
for j=1:numel(ecModel.rxns)
    % first get the relative transcript abundance
    if ~isempty(ecModel.rules{j})
        E(j) = splitGPRrule(ecModel.rules{j}, transcriptCounts);
    end
end
clear tmpTranscriptCounts


% adapt reaction bounds
ecModel.ub = E .* rxnKcats;

% keep upper bounds from model for exchange and transport reactions
transRxnIdx = ismember(model.rxns, findRxnsFromSubSystem(model,'Transport'));
noGprIdx = cellfun(@isempty, model.rules);
ecModel.ub(transRxnIdx|noGprIdx) = model.ub(transRxnIdx|noGprIdx);

% find NGAM reaction and keep upper bound 
ngamIdx = contains(model.rxnNames,'ngam');
if ~any(ngamIdx)
    error('no NGAM found using query ''ngam''')
end
ecModel.ub(ngamIdx)=model.ub(ngamIdx);

clear noGprIdx ngamIdx

% assign minimum upper limit to reactions with upper limit equal to zero
nzUbIdx = ecModel.ub>0;
ecModel.ub(~nzUbIdx) = min(ecModel.ub(nzUbIdx));

multFactor = 1E4;
fprintf('> all upper bounds were multiplied with %d! <\n',multFactor)
ecModel.ub = multFactor*ecModel.ub;

% report growth
s = optimizeCbModel(ecModel);
if ~isempty(s.x)&&s.f>0
    solution = s.x;
    growth = solution(ecModel.c==1);
    fprintf('\tobjective value: %.5f h^-1\n', growth)
else
    fprintf('\tno solution found!\n')
    growth = 0;
    solution = [];
end

end