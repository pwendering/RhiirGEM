% investigate growth response upon addition of myristate using the
% enzyme-constrained iRi1574 model
clear; clc
close all
% top-level directory
topDir = '';
% set solver for COBRA Toolbox
changeCobraSolver('ibm_cplex', 'all',0);

% load model
modelFile = [topDir, 'model/iRi1574.mat'];
load(modelFile)

%% remove unused genes and split reversible reactions
model = removeUnusedGenes(model);
model_irr = convertToIrreversible(model);
% additional carbon sources in the medium
cSources = {
    'r1533_e0' % D-Glucose
    'r1534_e0' % D-Fructose
    'r1006_e0' % myo-Inositol
    'r1002_e0' % Glycine
    };


%% kcats
% all kcats from BRENDA and SABIO-RK
maxKcatFile = fullfile(topDir,'data/kcats/kcat-reference-data.tsv');
% kcats associated to model reactions (will be created if it does not exist
% yet)
modelKcatsFile = fullfile(topDir, 'data/kcats/kcats_model.txt');
% associate maximum kcat per reaction
kcats = assign_kcats(model_irr, maxKcatFile, 'fungi');
% unknown kcats are assigned the median value of known values
kcats(kcats<=0) = median(kcats(kcats>0));
% convert to unit per hour
kcats = kcats * 3600; % [h^-1]

%% protein molecular weights
% UniProt results for R. irregularis
uniProtFile = fullfile(topDir,'data/uniprot.tab');
% calculate molecular weights for model proteins
[~, mw] = getAccountedProteinModel(model_irr, uniProtFile, []);
% convert from  g mol^-1 to g mmol^-1
mw = mw / 1000; % [g / mmol]

%% total protein content
C = 0.106; % g gDW^-1

%% Gurobi solver parameters
gurobiParams = struct;
gurobiParams.NumericFocus = 3;
gurobiParams.IntegralityFocus = 1;
gurobiParams.OutputFlag = 0;
gurobiParams.OptimalityTol = 1E-9;
gurobiParams.FeasibilityTol = 5e-7;
gurobiParams.IntFeasTol = 1e-9;

% define ratios between enzyme pools for myristate usage and the remaining
% reactions
a = linspace(0,1,1000);

% find reactions that a required for myristate usage
% as a proxy, use all peroxisomal enzymes
pxRxns = model_irr.rxns(~cellfun(@isempty,regexp(model_irr.rxns,'_x0[_fb]*$')));

% find associated genes
pxGenes = findGenesFromRxns(model_irr,pxRxns);
pxGenes = unique(vertcat(pxGenes{:}));
pxGeneIdx = findGeneIDs(model_irr,pxGenes);
remainGeneIdx = setdiff(1:numel(model_irr.genes),pxGeneIdx);

% convert problem to gurobi format
[~,problem] = enzymeFBA(model_irr,kcats,mw,C,[],false);
problem = cplex2gurobi(problem);

% find minimum uptake flux for palmitate at optimal growth
sol = gurobi(problem,gurobiParams);
bioMax = sol.objval;
tmpProblem = problem;
tmpProblem.lb(model_irr.c==1) = bioMax;
tmpProblem.obj(:) = 0;
tmpProblem.obj(findRxnIDs(model_irr, 'r1007_e0')) = 1;
tmpProblem.modelsense = 'min';
sol = gurobi(tmpProblem,gurobiParams);
minUptPalmitate = 0.1*sol.objval;
problem.ub(findRxnIDs(model_irr, 'r1007_e0')) = minUptPalmitate;

% do the same for additional carbon sources in the medium
cSourceIdx = cellfun(@(x)find(contains(model_irr.rxns,x)),cSources);
model_irr.lb(cSourceIdx) = 0;
tmpProblem.obj(:) = 0;
tmpProblem.obj(cSourceIdx) = 1;
tmpProblem.modelsense = 'min';
sol = gurobi(tmpProblem,gurobiParams);
problem.ub(cSourceIdx) = sol.x(cSourceIdx);
clear tmpProblem sol bioMax cSourceIdx

% remove original molecular crowding constraint
ccIdx = problem.rhs==C;
problem.A(ccIdx,:) = [];
problem.rhs(ccIdx,:) = [];
problem.sense(ccIdx) = [];

% loop over ratios, add the respective crowding constraints for the two
% enzyme pools
growthRatio = nan(size(a));
referenceGrowth = nan(size(a));
myrUptFlux = nan(size(a));
for i=1:numel(a)
    tmpProblem = problem;
    
    % divide protein pool
    pool_myr = a(i)*C;
    pool_remain = C-pool_myr;
    
    % add crowing constraints
    % myristate usage
    constr1 = zeros(1,size(tmpProblem.A,2));
    constr1(numel(model_irr.rxns)+pxGeneIdx) = mw(pxGeneIdx);
    tmpProblem.A = [tmpProblem.A; constr1];
    tmpProblem.rhs = [tmpProblem.rhs; pool_myr];
    tmpProblem.sense = [tmpProblem.sense, '<'];
    
    % remaining enzymes
    constr2 = zeros(1,size(tmpProblem.A,2));
    constr2(numel(model_irr.rxns)+remainGeneIdx) = mw(remainGeneIdx);
    tmpProblem.A = [tmpProblem.A; constr2];
    tmpProblem.rhs = [tmpProblem.rhs; pool_remain];
    tmpProblem.sense = [tmpProblem.sense, '<'];
    
    % as a reference, allow only for palmitate uptake
    tmpProblem.ub(findRxnIDs(model_irr,'r1633_e0')) = 0;
    s_ref = gurobi(tmpProblem,gurobiParams);
    
    % solve with addition of myristate
    tmpProblem.ub(findRxnIDs(model_irr,'r1633_e0')) = 1000;
    s = gurobi(tmpProblem,gurobiParams);
    
    if isequal(s_ref.status,'OPTIMAL') && isequal(s.status,'OPTIMAL')
        if s_ref.objval>0 && s.objval>0
            fprintf('Myristate uptake: %.5g\n', s.x(findRxnIDs(model_irr,'r1633_e0')))
            referenceGrowth(i) = s_ref.objval;
            growthRatio(i) = s.objval/s_ref.objval;
            myrUptFlux(i) = s.x(findRxnIDs(model_irr,'r1633_e0'));
        end
    end

end
% smooth_growth_ratio = smooth(growthRatio,50);
% smooth_reference_growth = smooth(referenceGrowth,50);
yyaxis left
plot(a,referenceGrowth,'LineWidth', 1.5)
% plot(a,smooth_reference_growth,'LineWidth', 1.5)
hold on
% plot(a,myrUptFlux,'LineWidth', 1.5)
ylabel('\mu_{ref} [h^{-1}]')
xlabel('Ratio between peroxisomal and remaining enzyme pool')

yyaxis right
plot(a,growthRatio,'LineWidth', 1.5)
% plot(a,smooth_growth_ratio,'LineWidth', 1.5)

ylabel('Z^{\alpha}')

box on
set(gca,'LineWidth', 1.3, 'FontSize', 14)

print([topDir 'results/figures/Figure-S4.png'], '-painters', '-dpng')

%% check ratio in expression data
transcriptomicDir = fullfile(topDir, 'data/transcriptomic-data');
transcriptomicFiles = {...
    fullfile(transcriptomicDir, 'transcriptomic_ERM.csv'),...
    fullfile(transcriptomicDir, 'transcriptomic_HYP.csv'),...
    fullfile(transcriptomicDir, 'transcriptomic_ARB.csv'),...
    };

for i=1:numel(transcriptomicFiles)
    transcriptCounts = readtable(transcriptomicFiles{i});
    
    % take the subset of genes in the model
    idx = cell2mat(cellfun(@(x)find(ismember(transcriptCounts.Var1, x)), model.genes, 'un', 0));
    transcriptCounts = transcriptCounts.Var2(idx);
    % scale counts to get relative abundances
    transcriptCounts = transcriptCounts / sum(transcriptCounts);
    
    pxFraction = sum(transcriptCounts(pxGeneIdx));
    fprintf('Peroxisome fraction: %.2g\n', pxFraction)
       
end